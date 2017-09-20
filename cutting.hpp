#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_CUTTING_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_CUTTING_HPP 1

// ell.cpp
#include <cmath>
#include <tuple>
#include <valarray>

// ell = { x | (x - xc)' * P^-1 * (x - xc) <= 1 }
class ell {
private:
  using Vec = std::valarray<double>;
  // using Mat = std::valarray<Vec>;

  class Mat : public std::valarray<Vec> {
  public:
    Mat(const Vec& v, size_t n) : std::valarray<Vec>(v, n) {}
    Vec operator*(const Vec &g) const {
      size_t n = g.size();
      Vec gt(n);
      for (size_t i = 0; i < n; ++i) {
        gt[i] = ((*this)[i] * g).sum();
      }
      return gt;
    }
  };

  size_t n;
  Mat _P;
  Vec _xc;
  double _c1;

protected:
  void update_P(double sig, double delta, const Vec &gt) {
    for (size_t i = 0; i < n; ++i) {
      double temp = sig * gt[i];
      for (size_t j = i; j < n; ++j) {
        _P[i][j] -= temp * gt[j];
        _P[i][j] *= delta;
      }
    }

    // Make symmetric
    for (size_t i = 0; i < n - 1; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        _P[j][i] = _P[i][j];
      }
    }
  }

public:
  ell(double val, const Vec &x)
      : n{x.size()},
        _P(Mat(Vec(0.0, n), n)), _xc{x}, _c1{n * n / (n * n - 1.0)} {
    for (size_t i = 0; i < n; ++i) {
      _P[i][i] = val;
    }
  }

  ell(const Vec &r, const Vec &x)
      : n{x.size()},
        _P(Mat(Vec(0.0, n), n)), _xc{x}, _c1{n * n / (n * n - 1.0)} {
    for (size_t i = 0; i < n; ++i) {
      _P[i][i] = r[i];
    }
  }

  Vec &xc() { return _xc; }

  double update(const Vec &g) { //  central cut
    auto Pg = _P * g;
    double tsq = (g * Pg).sum();
    double rho = 1.0 / (n + 1);
    double tau = sqrt(tsq);
    _xc -= (rho / tau) * Pg;
    update_P(2 * rho / tsq, _c1, Pg);
    return tau;
  }

  std::tuple<int, double> update(const Vec &g, double beta) { //  deep cut
    auto Pg = _P * g;
    double tsq = (g * Pg).sum();
    double tau = sqrt(tsq);

    if (beta > tau) { // no sol'n
      return {1, tau};
    }
    if (n * beta < -tau) { // no effect
      return {3, tau};
    }

    double alpha = beta / tau;
    double rho = (1 + n * alpha) / (n + 1);
    double sigma = 2 * rho / (1 + alpha);
    double delta = _c1 * (1 - alpha * alpha);
    _xc -= (rho / tau) * Pg;
    update_P(sigma / tsq, delta, Pg);
    return {0, tau};
  }

  std::tuple<int, double> update(const Vec &g,
                                 const Vec &beta) { // parallel or deep cut
    auto Pg = _P * g;
    double tsq = (g * Pg).sum();
    double tau = sqrt(tsq);
    auto alpha = beta / tau;

    if (alpha[0] > 1) { // no sol'n
      return {1, tau};
    }

    double rho, sigma, delta;

    if (alpha.size() == 1 || alpha[1] >= 1.0) { // deep cut
      double a0 = alpha[0];
      if (n * a0 < -1.0) { // no effect
        return {3, tau};
      }
      rho = (1 + n * a0) / (n + 1);
      sigma = 2 * rho / (1 + a0);
      delta = _c1 * (1 - a0 * a0);
    } else { // parallel cut
      double aprod = alpha[0] * alpha[1];
      if (n * aprod < -1.0) { // no effect
        return std::make_tuple(3, tau);
      }
      auto asq = alpha * alpha;
      double asum = alpha[0] + alpha[1];
      double asqdiff = asq[1] - asq[0];
      double xi =
          sqrt(4 * (1 - asq[0]) * (1 - asq[1]) + (n * n) * asqdiff * asqdiff);
      sigma = (n + (2 * (1 + aprod - xi / 2) / (asum * asum))) / (n + 1);
      rho = asum * sigma / 2;
      delta = _c1 * (1 - (asq[0] + asq[1] - xi / n) / 2);
    }
    _xc = _xc - (rho / tau) * Pg;
    update_P(sigma / tsq, delta, Pg);
    return {0, tau};
  }
}; // } ell

// -- Generic Cutting-plane method for solving convex optimization problem
//
// input
//         oracle        perform assessment on x0
//          E(P,xc)       ellipsoid containing x*
//         t             best-so-far optimal sol'n
//         max_it        maximum number of iterations
//         tol           error tolerance
//
// output
//         x             solution vector
//         iter          number of iterations performed
template <class F>
auto cutting_plane_dc(F &access, ell &E, double t, int max_it, double tol) {
  using Vec = std::valarray<double>;
  Vec x_best = E.xc();
  int iter = 1, flag = 0;

  for (; iter < max_it; ++iter) {
    auto [g, h, t1] = access(E.xc(), t);
    if (t != t1) { // best t obtained
      flag = 1;
      t = t1;
      x_best = E.xc();
    }
    auto [status, tau] = E.update(g, h);
    if (status == 1) {
      break;
    }
    if (tau < tol) {
      status = 2;
      break;
    } // no more,
  }
  return std::make_tuple(x_best, t, iter, flag, 0);
} // END

template <class F>
auto cutting_plane_q(F &access, ell &E, double t, int max_it, double tol) {
  using Vec = std::valarray<double>;
  auto norm = [](const Vec &v) -> double { return sqrt((v * v).sum()); };
  Vec x_last = E.xc();
  Vec x_best = E.xc();
  int iter = 1, flag = 0;
  int count = 20;

  for (; iter < max_it; ++iter) {
    auto [g, h, t1, x, loop] = access(E.xc(), t);
    if (loop == 1) {
      h += (g * (x - E.xc())).sum();
    }

    if (t != t1) { // best t obtained
      flag = 1;
      t = t1;
      x_best = E.xc();
    }
    auto [status, tau] = E.update(g, h);

    if (status == 1) {
      break;
    } 
    if (status == 3) { // retry 20 times
      --count;
      if (count == 0) {
        if (flag == 0) {
          x_best = x;
        }
        break;
      }
      continue;
    }
    
    if (tau < tol) {
      status = 2;
      break;
    } // no more,
    if (norm(x_last - E.xc()) < 1e-8) {
      status = 4;
      break;
    }

    x_last = E.xc();
    count = 20; // restart the count
  }

  return std::make_tuple(x_best, t, iter, flag, 0);
} // END

#endif