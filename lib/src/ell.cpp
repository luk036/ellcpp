#include <ell.hpp>
#include <cmath>
#include <tuple>

/**
 * @brief Central Cut
 */
ell::params_t ell::calc_cc(std::size_t n) const {
  double rho = 1. / (n + 1);
  double sigma = 2. * rho;
  double delta = _c1;
  return std::tuple{rho, sigma, delta};
}

/**
 * @brief Deep Cut
 */
ell::return_t ell::calc_dc(double alpha) const {
  std::size_t n = _xc.size();
  if (alpha == 0.) {
    return std::tuple{0, this->calc_cc(n)};
  }
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

/** Situation when feasible cut. */
ell::params_t ell::calc_ll_cc(double a1, std::size_t n) const {
  double asq1 = a1 * a1;
  double nasq1 = n * asq1;
  double xi = std::sqrt(nasq1 * nasq1 - 4. * asq1 + 4.);
  double sigma = (n + (2. - xi) / asq1) / (n + 1);
  double rho = a1 * sigma / 2.;
  double delta = _c1 * (1 - (asq1 - xi / n) / 2.);
  return std::tuple{rho, sigma, delta};
}

ell::params_t ell::calc_ll_general(double a0, double a1, std::size_t n) const {
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


ell1d::return_t ell1d::update(double g, double beta) {
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
