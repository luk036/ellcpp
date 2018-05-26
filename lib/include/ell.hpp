#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ELL_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ELL_HPP 1

#include <cmath>
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

public:
  bool _use_parallel = true;

private:
  std::size_t n;
  double _c1;
  double _kappa;
  Vec _xc;
  Mat _Q;

public:
  template <typename T, typename V>
  ell(const T &val, const V &x)
      : n{x.size()}, _c1{n * n / (n * n - 1.0)}, _xc{x} {
    if constexpr (std::is_scalar<T>::value) { // C++17
      _Q = xt::eye(n);
      _kappa = val;
    } else {
      _Q = xt::diag(val);
      _kappa = 1.;
    }
  }

  auto &xc() { return _xc; }

  ell(const ell &E) = default;

  template <typename T, typename V>
  auto update_core(const V &g, const T &beta) {
    Vec Qg = xt::linalg::dot(_Q, g);
    auto tsq = xt::linalg::dot(g, Qg)();
    auto tau = std::sqrt(_kappa * tsq);
    auto alpha = beta / tau;
    auto [status, rho, sigma, delta] = this->calc_ll(alpha);
    if (status != 0) {
      return std::tuple{status, tau}; // g++-7 is ok with this
    }
    _xc -= (_kappa * rho / tau) * Qg;
    _Q -= (sigma / tsq) * xt::linalg::outer(Qg, Qg);
    // _Q *= delta;
    _kappa *= delta;
    return std::tuple{status, tau}; // g++-7 is ok
    // return std::pair{status, tau}; // workaround for clang++ 6
  }

  /**
   * @brief Central Cut
   */
  auto calc_cc() {
    auto n = _xc.size();
    auto rho = 1. / (n + 1);
    auto sigma = 2. * rho;
    auto delta = _c1;
    return std::tuple{0, rho, sigma, delta};
  }

  /**
   * @brief Deep Cut
   */
  auto calc_dc(double alpha) {
    if (alpha == 0.) {
      return this->calc_cc();
    }
    auto n = _xc.size();
    // auto [status, rho, sigma, delta] = std::tuple{0, 0., 0., 0.0};
    auto status = 0;
    auto rho = 0., sigma = 0., delta = 0.;

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
      auto a0 = alpha[0];
      if (alpha.shape()[0] < 2 || !_use_parallel) {
        return this->calc_dc(a0);
      }
      auto a1 = alpha[1];
      // auto [a0, a1] = alpha;
      if (a1 >= 1.) {
        return this->calc_dc(a0);
      }
      auto n = _xc.size();

      // auto [status, rho, sigma, delta] = std::tuple{0, 0.0, 0.0, 0.0};
      auto status = 0;
      auto rho = 0., sigma = 0., delta = 0.;
      auto aprod = a0 * a1;
      if (a0 > a1) {
        status = 1; // no sol'n
      } else if (n * aprod < -1.) {
        status = 3; // no effect
      } else {
        // auto asq = bnu::element_prod(alpha, alpha);
        auto asq = alpha * alpha;
        auto asum = a0 + a1;
        // auto [asq0, asq1] = asq;
        auto asq0 = asq[0], asq1 = asq[1];
        auto asqdiff = asq1 - asq0;
        auto xi = std::sqrt(4. * (1. - asq0) * (1. - asq1) +
                            n * n * asqdiff * asqdiff);
        sigma = (n + (2. * (1. + aprod - xi / 2.) / (asum * asum))) / (n + 1);
        rho = asum * sigma / 2.;
        delta = _c1 * (1. - (asq0 + asq1 - xi / n) / 2.);
      }
      return std::tuple{status, rho, sigma, delta};
    }
  }

  template <typename T> auto update(const Vec &g, const T &beta) {
    return this->update_core(g, beta);
  }

}; // } ell

class ell1d {
private:
  double _r;
  double _xc;

public:
  ell1d(double l, double u) : _r{(u - l) / 2}, _xc{l + _r} {}

  auto &xc() { return _xc; }

  ell1d(const ell1d &E) = default;

  auto update(double g, double beta) {
    auto tau = std::abs(_r * g);
    if (beta == 0.) {
      _r /= 2;
      if (g > 0.) {
        _xc -= _r;
      } else {
        _xc += _r;
      }
      return std::tuple{0, tau};
    }
    if (beta > tau) {
      return std::tuple{1, tau}; // no sol'n
    }
    if (beta < -tau) {
      return std::tuple{3, tau}; // no effect
    }
    double l, u;
    double bound = _xc - beta / g;
    if (g > 0.) {
      u = bound;
      l = _xc - _r;
    } else {
      l = bound;
      u = _xc + _r;
    }
    _r = (u - l) / 2;
    _xc = l + _r;
    return std::tuple{0, tau};
  }

}; // } ell1d

#endif