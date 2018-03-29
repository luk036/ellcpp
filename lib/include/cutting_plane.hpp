#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_CUTTING_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_CUTTING_PLANE_HPP 1

#include <cmath>
#include <tuple>

/**
 * @brief Cutting-plane method for solving convex feasibility problem
 * 
 * @tparam Oracle 
 * @tparam Space 
 * @tparam T 
 * @param[in] assess   perform assessment on x0
 * @param[in] S        search Space containing x*
 * @param[in] t        best-so-far optimal sol'n
 * @param[in] max_it   maximum number of iterations
 * @param[in] tol      error tolerance
 * @return x      solution vector 
 * @return iter   number of iterations performed
 * @return flag   solution found or not 
 * @return status how is the final cut 
 */
template <typename Oracle, typename Space, typename T>
auto cutting_plane_feas(Oracle &assess, Space &S, T t, int max_it = 1000,
                        double tol = 1e-8) {
  auto flag = 0;
  auto iter = 1, status = 0;
  for (; iter <= max_it; ++iter) {
    auto [g, h, t1] = assess(S.xc(), t);
    if (t != t1) {
      flag = 1;
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
  return std::tuple{S.xc(), iter, flag, status};
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
**/
template <typename Oracle, typename Space, typename T>
auto cutting_plane_dc(Oracle &assess, Space &S, T t, int max_it = 1000,
                      double tol = 1e-8) {
  auto x_best = S.xc();
  auto flag = 0;
  auto iter = 1, status = 0;
  for (; iter <= max_it; ++iter) {
    auto [g, h, t1] = assess(S.xc(), t);
    if (t != t1) { // best t obtained
      flag = 1;
      t = t1;
      x_best = S.xc();
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
  return std::tuple{x_best, t, iter, flag, status};
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
// #include <boost/numeric/ublas/symmetric.hpp>
// namespace bnu = boost::numeric::ublas;
#include <xtensor/xarray.hpp>
#include <xtensor-blas/xlinalg.hpp>

template <typename Oracle, typename Space, typename T>
auto cutting_plane_q(Oracle &assess, Space &S, T t, int max_it = 1000,
                     double tol = 1e-8) {
  auto flag = 0;
  auto x_best = S.xc();
  auto iter = 1, status = 1;
  for (; iter < max_it; ++iter) {
    auto [g, h, t1, x, loop] = assess(S.xc(), t, (status != 3) ? 0 : 1);
    if (status != 3) {
      if (loop == 1) { // discrete sol'n
        h += xt::linalg::dot(g, x - S.xc())();
      }
    } else { // can't cut in the previous iteration
      if (loop == 0) {
        break; // no more alternative cut
      }
      h += xt::linalg::dot(g, x - S.xc())();
    }
    if (t != t1) { // best t obtained
      flag = 1;
      t = t1;
      x_best = x;
    }
    double tau;
    std::tie(status, tau) = S.update(g, h);
    if (status == 1) {
      break;
    }
    if (tau < tol) {
      status = 2;
      break;
    }
  }

  return std::tuple{x_best, t, iter, flag, status};
} // END

#endif