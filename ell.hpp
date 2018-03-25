#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ELL_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ELL_HPP 1

#include <tuple>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

/**
 * @brief Ellipsoid Search Space
 *
 * ell = { x | (x - xc)' * P^-1 * (x - xc) <= 1 }
 */
class ell {
  using Mat = xt::xarray<double>;
  using Vec = xt::xarray<double>;

private:
  size_t n;
  double _c1;
  Vec _xc;
  Mat _P;

public:
  template <typename T, typename V>
  ell(const T &val, const V &x)
      : n{x.size()}, _c1{n * n / (n * n - 1.0)}, _xc{x} {
    if constexpr (std::is_scalar<T>::value) { // C++17
      _P = val * xt::eye(n);
    } else {
      _P = xt::diag(val);
    }
  }

  auto &xc() { return _xc; }

  template <typename T, typename V>
  auto update_core(const V &g, const T &beta) {
    Vec Pg = xt::linalg::dot(_P, g);
    auto tsq = xt::linalg::dot(g, Pg)();
    auto tau = std::sqrt(tsq);
    auto alpha = beta / tau;
    auto [status, rho, sigma, delta] = this->calc_ll(alpha);
    if (status != 0) {
      return std::tuple{status, tau}; // g++-7 is ok with this
    }
    _xc -= (rho / tau) * Pg;
    _P -= (sigma / tsq) * xt::linalg::outer(Pg, Pg);
    _P *= delta;
    return std::tuple{status, tau}; // g++-7 is ok
    // return std::pair{status, tau}; // workaround for clang++ 6
  }

  /**
   * @brief Central Cut
   */
  auto calc_cc() {
    auto n = _xc.size();
    auto rho = 1.0 / (n + 1);
    auto sigma = 2.0 * rho;
    auto delta = _c1;
    return std::tuple{0, rho, sigma, delta};
  }

  /**
   * @brief Deep Cut
   */
  auto calc_dc(double alpha) {
    if (alpha == 0.0) {
      return this->calc_cc();
    }
    auto n = _xc.size();
    auto [status, rho, sigma, delta] = std::tuple{0, 0.0, 0.0, 0.0};
    if (alpha > 1.) {
      status = 1; // no sol'n
    } else if (n * alpha < -1.) {
      status = 3; // no effect
    } else {
      rho = (1.0 + n * alpha) / (n + 1);
      sigma = 2.0 * rho / (1.0 + alpha);
      delta = _c1 * (1.0 - alpha * alpha);
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
      auto n = _xc.size();

      // auto [status, rho, sigma, delta] = std::tuple{0, 0.0, 0.0, 0.0};
      auto status = 0;
      auto rho = 0.0, sigma = 0.0, delta = 0.0;
      auto aprod = a0 * a1;
      if (a0 > a1) {
        status = 1; // no sol'n
      } else if (n * aprod < -1.0) {
        status = 3; // no effect
      } else {
        // auto asq = bnu::element_prod(alpha, alpha);
        auto asq = alpha * alpha;
        auto asum = a0 + a1;
        auto asqdiff = asq(1) - asq(0);
        auto xi = std::sqrt(4.0 * (1.0 - asq(0)) * (1.0 - asq(1)) +
                            n * n * asqdiff * asqdiff);
        sigma =
            (n + (2.0 * (1.0 + aprod - xi / 2.0) / (asum * asum))) / (n + 1);
        rho = asum * sigma / 2.0;
        delta = _c1 * (1.0 - (asq(0) + asq(1) - xi / n) / 2.0);
      }
      return std::tuple{status, rho, sigma, delta};
    }
  }

  template <typename T> auto update(const Vec &g, const T &beta) {
    return this->update_core(g, beta);
  }

}; // } ell

#endif