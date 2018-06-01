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
  using Arr = xt::xarray<double>;

public:
  bool _use_parallel = true;

private:
  std::size_t n;
  double _c1;
  double _kappa;
  Arr _xc;
  Arr _Q;

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

  void set_xc(const Arr &xc) { _xc = xc; }

  ell(const ell &E) = default;

  template <typename T, typename V>
  auto update_core(const V &g, const T &beta) {
    Arr Qg = xt::linalg::dot(_Q, g);
    auto tsq = xt::linalg::dot(g, Qg)();
    auto tau = std::sqrt(_kappa * tsq);
    auto alpha = beta / tau;
    auto [status, params] = this->calc_ll(alpha);
    if (status != 0) {
      return std::tuple{status, tau}; // g++-7 is ok with this
    }
    auto [rho, sigma, delta] = params;
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
    return std::tuple{rho, sigma, delta};
  }

  /**
   * @brief Deep Cut
   */
  auto calc_dc(double alpha) {
    if (alpha == 0.) {
      return std::tuple{0, this->calc_cc()};
    }
    auto n = _xc.size();
    // auto [status, rho, sigma, delta] = std::tuple{0, 0., 0., 0.0};
    auto status = 0;
    auto params = std::tuple{0., 0., 0.};
    
    if (alpha > 1.) {
      status = 1; // no sol'n
    } else if (n * alpha < -1.) {
      status = 3; // no effect
    } else {
      double rho = (1. + n * alpha) / (n + 1);
      double sigma = 2. * rho / (1. + alpha);
      double delta = _c1 * (1. - alpha * alpha);
      params = std::tuple{rho, sigma, delta};
    }
    return std::tuple{status, params};
  }

  /* parallel or deep cut */
  template <typename T> auto calc_ll(const T &alpha) {
    if constexpr (std::is_scalar<T>::value) { // C++17
      return this->calc_dc(alpha);
    } else { // parallel cut
      auto a0 = alpha[0];
      if (alpha.shape()[0] < 2 || !_use_parallel) {
        return this->calc_dc(a0);
      }
      auto a1 = alpha[1];
      if (a1 >= 1.) {
        return this->calc_dc(a0);
      }
      auto n = _xc.size();
      auto status = 0;
      auto params = std::tuple{0., 0., 0.};
      if (a0 > a1) {
        status = 1; // no sol'n
      } else if (n * a0 * a1 < -1.) {
        status = 3; // no effect
      } else if (a0 == 0.) {
        params = this->calc_ll_cc(a1, n);
      } else {
        params = this->calc_ll_general(a0, a1, n);
      }
      return std::tuple{status, params};
    }
  }

  /** Situation when feasible cut. */
  auto calc_ll_cc(double a1, std::size_t n) {
    double asq1 = a1 * a1;
    double nasq1 = n * asq1;
    double xi = std::sqrt(nasq1 * nasq1 - 4. * asq1 + 4.);
    double sigma = (n + (2. - xi) / asq1) / (n + 1);
    double rho = a1 * sigma / 2.;
    double delta = _c1 * (1 - (asq1 - xi / n) / 2.);
    return std::tuple{rho, sigma, delta};
  }

  auto calc_ll_general(double a0, double a1, std::size_t n) {
    double asum = a0 + a1;
    double asq0 = a0 * a0, asq1 = a1 * a1;
    double asqdiff = asq1 - asq0;
    double nasqdiff = n * asqdiff;
    double xi = std::sqrt(4. * (1. - asq0) * (1. - asq1) + nasqdiff * nasqdiff);
    double sigma =
        (n + (2. * (1. + a0 * a1 - xi / 2.) / (asum * asum))) / (n + 1);
    double rho = asum * sigma / 2.;
    double delta = _c1 * (1. - (asq0 + asq1 - xi / n) / 2.);
    return std::tuple{rho, sigma, delta};
  }

  template <typename T> auto update(const Arr &g, const T &beta) {
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