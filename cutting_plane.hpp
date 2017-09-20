#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_CUTTING_PLANE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_CUTTING_PLANE_HPP 1

#include <cmath>
#include <tuple>

/**
    Cutting-plane method for solving convex feasibility problem
    input
             assess        perform assessment on x0
             S(xc)         Search Space containing x*
             t             best-so-far optimal sol'n
             max_it        maximum number of iterations
             tol           error tolerance
    output
             x             solution vector
             iter          number of iterations performed
**/
template <typename T1, typename T2, typename T3>
auto cutting_plane_feas(T1 &assess, T2 &S, T3 t, int max_it = 1000,
                        double tol = 1e-8) {
  int iter = 1, flag = 0, status = 0;
  for (; iter <= max_it; ++iter) {
    auto [g, h, flag] = assess(S.xc(), t);
    if (flag) {
      break;
    }
    double tau;
    std::tie(status, tau) = S.update(g, h);
    if (status == 1) {
      break;
    }
    if (tau < tol) { // no more
      status = 2;
      break;
    }
  }
  return std::make_tuple(S.xc(), iter, flag, status);
}

/**
    Cutting-plane method for solving convex optimization problem
    input
             assess        perform assessment on x0
             S(xc)         Search Space containing x*
             t             initial best-so-far value
             max_it        maximum number of iterations
             tol           error tolerance
    output
             x_best        solution vector
             t             best-so-far optimal value
             iter          number of iterations performed
    '''
**/
template <typename T1, typename T2, typename T3>
auto cutting_plane_dc(T1 &assess, T2 &S, T3 t, int max_it = 1000,
                      double tol = 1e-8) {
  auto x_best = S.xc();
  int iter = 1, flag = 0, status = 0;
  double tau;
  for (; iter <= max_it; ++iter) {
    auto [g, h, t1] = assess(S.xc(), t);
    if (t != t1) { // best t obtained
      flag = 1;
      t = t1;
      x_best = S.xc();
    }
    std::tie(status, tau) = S.update(g, h);
    if (status == 1) {
      break;
    }
    if (tau < tol) { // no more
      status = 2;
      break;
    }
  }
  return std::make_tuple(x_best, t, iter, flag, status);
} // END

/**
    Cutting-plane method for solving convex discrete optimization problem
    input
             oracle        perform assessment on x0
             S(xc)         Search space containing x*
             t             best-so-far optimal sol'n
             max_it        maximum number of iterations
             tol           error tolerance
    output
             x             solution vector
             iter          number of iterations performed
**/
#include <boost/numeric/ublas/symmetric.hpp>
namespace bnu = boost::numeric::ublas;

template <typename T1, typename T2, typename T3>
auto cutting_plane_q(T1 &assess, T2 &S, T3 t, int max_it = 1000,
                     double tol = 1e-8) {
  auto x_best = S.xc();
  int iter = 1, flag = 0, status = 1;
  double tau;

  for (; iter < max_it; ++iter) {
    auto[g, h, t1, x, loop] = assess(S.xc(), t, (status != 3) ? 0 : 1);
    if (status != 3) {
      if (loop == 1) { // discrete sol'n
        h += bnu::inner_prod(g, x - S.xc());
      }
    } else { // can't cut in the previous iteration
      if (loop == 0) {
        break; // no more alternative cut
      }
      h += bnu::inner_prod(g, x - S.xc());
    }
    if (t != t1) { // best t obtained
      flag = 1;
      t = t1;
      x_best = x;
    }
    std::tie(status, tau) = S.update(g, h);
    if (status == 1) {
      break;
    }
    if (tau < tol) {
      status = 2;
      break;
    }
  }

  return std::make_tuple(x_best, t, iter, flag, status);
} // END

#endif