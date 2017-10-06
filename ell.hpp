#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_ELL_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_ELL_HPP 1

#include <boost/numeric/ublas/symmetric.hpp>
#include <tuple>

namespace bnu = boost::numeric::ublas;

// ell = { x | (x - xc)' * P^-1 * (x - xc) <= 1 }
class ell {
  using Mat = bnu::symmetric_matrix<double, bnu::upper>;
  using Vec = bnu::vector<double>;

private:
  size_t n;
  Mat _P;
  Vec _xc;
  double _c1;

public:
  template <typename T, typename V>
  ell(const T &val, const V &x)
      : n{x.size()}, _P(n, n), _xc{x}, _c1{n * n / (n * n - 1.0)} {
    for (auto i = 0U; i < n; ++i) {
      if constexpr (std::is_scalar<T>::value) { // C++17
        _P(i, i) = val;
      } else {
        _P(i, i) = val(i);
      }
    }
  }

  auto &xc() { return _xc; }

  template <typename T, typename V>
  auto update_core(const V &g, const T &beta) {
    Vec Pg = bnu::prod(this->_P, g);
    auto tsq = bnu::inner_prod(g, Pg);
    double tau = std::sqrt(tsq);
    auto alpha = beta / tau;
    auto [status, rho, sigma, delta] = this->calc_ll(alpha);
    if (status == 0) {
      this->_xc -= (rho / tau) * Pg;
      this->_P -= (sigma / tsq) * bnu::outer_prod(Pg, Pg);
      this->_P *= delta;
    }
    //return std::tuple{status, tau}; // g++-7 is ok with this
    return std::pair{status, tau}; // workaround for clang++ 6
  }

  auto calc_cc() {
    /* central cut */
    auto n = this->_xc.size();
    auto rho = 1.0 / (n + 1);
    auto sigma = 2.0 * rho;
    auto delta = this->_c1;
    return std::tuple{0, rho, sigma, delta};
  }

  auto calc_dc(double alpha) {
    /* deep cut */
    if (alpha == 0.0) {
      return this->calc_cc();
    }
    auto n = this->_xc.size();
    auto [status, rho, sigma, delta] = std::tuple{0, 0.0, 0.0, 0.0};
    if (alpha > 1.) {
      status = 1; // no sol'n
    } else if (n * alpha < -1.) {
      status = 3; // no effect
    } else {
      rho = (1.0 + n * alpha) / (n + 1);
      sigma = 2.0 * rho / (1.0 + alpha);
      delta = this->_c1 * (1.0 - alpha * alpha);
    }
    return std::tuple{status, rho, sigma, delta};
  }

  template <typename T> auto calc_ll(const T &alpha) {
    /* parallel or deep cut */
    if constexpr (std::is_scalar<T>::value) { // C++17
      return this->calc_dc(alpha);
    } else { // parallel cut
      // auto a0 = alpha(0), a1 = alpha(1);
      auto [a0, a1] = alpha;
      if (a1 >= 1.0) {
        return this->calc_dc(a0);
      }
      auto n = this->_xc.size();
      auto [status, rho, sigma, delta] = std::tuple{0, 0.0, 0.0, 0.0};
      auto aprod = a0 * a1;
      if (a0 > a1) {
        status = 1; // no sol'n
      } else if (n * aprod < -1.0) {
        status = 3; // no effect
      } else {
        auto asq = bnu::element_prod(alpha, alpha);
        auto asum = a0 + a1;
        auto asqdiff = asq(1) - asq(0);
        auto xi = std::sqrt(4.0 * (1.0 - asq(0)) * (1.0 - asq(1)) +
                            n * n * asqdiff * asqdiff);
        sigma =
            (n + (2.0 * (1.0 + aprod - xi / 2.0) / (asum * asum))) / (n + 1);
        rho = asum * sigma / 2.0;
        delta = this->_c1 * (1.0 - (asq(0) + asq(1) - xi / n) / 2.0);
      }
      return std::tuple{status, rho, sigma, delta};
    }
  }

  template <typename T> auto update(const Vec &g, const T &beta) {
    return this->update_core(g, beta);
  }

}; // } ell

#endif