#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_CUTTING_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_CUTTING_PLANE_HPP 1

#include <cmath>
#include <tuple>
// #include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

/**
 * @brief Options
 *
 */
struct Options {
    unsigned int max_it = 2000;
    double tol = 1e-8;
};

/**
 * @brief
 *
 * @tparam Oracle
 * @tparam Space
 * @param evaluate
 * @param I
 * @param options
 * @return auto
 */
template <typename Oracle, typename Space>
auto bsearch(Oracle &evaluate, Space &I, const Options &options = Options()) {
    // assume monotone
    bool feasible = false;
    auto &[l, u] = I;
    auto t = l + (u - l) / 2;
    auto niter = 1U;
    for (; niter <= options.max_it; ++niter) {
        if (evaluate(t)) { // feasible sol'n obtained
            feasible = true;
            u = t;
        } else {
            l = t;
        }
        auto tau = (u - l) / 2;
        t = l + tau;
        if (tau < options.tol) {
            break;
        }
    }
    return std::tuple{u, niter, feasible};
}

/**
 * @brief
 *
 * @tparam Oracle
 * @tparam Space
 */
template <typename Oracle, typename Space> //
class bsearch_adaptor {
  private:
    Oracle &_P;
    Space &_S;
    Options _options;

  public:
    explicit bsearch_adaptor(Oracle &P, Space &S,
                             const Options &options = Options())
        : _P{P}, //
          _S{S}, //
          _options{options} {}

    auto x_best() const { return _S.xc(); }

    auto operator()(double t) {
        Space S(_S);
        _P.update(t);
        auto [x, _1, feasible, _2] = cutting_plane_feas(_P, S, _options);
        if (feasible) {
            _S.xc() = x;
            return true;
        }
        return false;
    }
};

/**
 * @brief Cutting-plane method for solving convex feasibility problem
 *
 * @tparam Oracle
 * @tparam Space
 * @tparam T
 * @param[in] evaluate   perform assessment on x0
 * @param[in] S        search Space containing x*
 * @param[in] max_it   maximum number of iterations
 * @param[in] tol      error tolerance
 * @return x      solution vector
 * @return niter  number of iterations performed
 * @return feasible   solution found or not
 * @return status how is the final cut
 */
template <typename Oracle, typename Space>
auto cutting_plane_feas(Oracle &evaluate, Space &S,
                        const Options &options = Options()) {
    bool feasible = false;
    auto niter = 1U, status = 0U;

    for (; niter <= options.max_it; ++niter) {
        auto [g, h, flag] = evaluate(S.xc());
        if (flag) { // feasible sol'n obtained
            feasible = true;
            break;
        }
        double tsq;
        std::tie(status, tsq) = S.update(g, h);
        if (status != 0) {
            break;
        }
        if (tsq < options.tol) { // no more
            status = 2;
            break;
        }
    }
    return std::tuple{S.xc(), niter, feasible, status};
}

/**
 * @brief Cutting-plane method for solving convex problem
 *
 * @tparam Oracle
 * @tparam Space
 * @tparam T
 * @param[in] evaluate   perform assessment on x0
 * @param[in] S        search Space containing x*
 * @param[in] t        best-so-far optimal sol'n
 * @param[in] max_it   maximum number of iterations
 * @param[in] tol      error tolerance
 * @return x      solution vector
 * @return niter   number of iterations performed
 * @return feasible   solution found or not
 * @return status how is the final cut
 */
template <typename Oracle, typename Space, typename T>
auto cutting_plane_dc(Oracle &evaluate, Space &S, T t,
                      const Options &options = Options()) {
    using Arr = xt::xarray<double>;

    bool feasible = false;
    auto x_best = Arr{S.xc()};
    auto niter = 1U, status = 0U;
    for (; niter <= options.max_it; ++niter) {
        auto [g, h, t1] = evaluate(S.xc(), t);
        if (t != t1) { // best t obtained
            feasible = true;
            t = t1;
            x_best = S.xc();
        }
        double tsq;
        std::tie(status, tsq) = S.update(g, h);
        if (status == 1) {
            break;
        }
        if (tsq < options.tol) { // no more
            status = 2;
            break;
        }
    }
    return std::tuple{std::move(x_best), t, niter, feasible, status};
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
             niter          number of iterations performed
**/
// #include <boost/numeric/ublas/symmetric.hpp>
// namespace bnu = boost::numeric::ublas;
#include <xtensor-blas/xlinalg.hpp>
// #include <xtensor/xarray.hpp>

/**
 * @brief
 *
 * @tparam Oracle
 * @tparam Space
 * @tparam T
 * @param evaluate
 * @param S
 * @param t
 * @param options
 * @return auto
 */
template <typename Oracle, typename Space, typename T>
auto cutting_plane_q(Oracle &evaluate, Space &S, T t,
                     const Options &options = Options()) {
    using Arr = xt::xarray<double>;

    bool feasible = false;
    auto x_best = Arr{S.xc()};
    auto status = 1U;
    auto niter = 1U;
    for (; niter < options.max_it; ++niter) {
        auto [g, h, t1, x0, loop] = evaluate(S.xc(), t, (status != 3) ? 0 : 1);
        if (status != 3) {
            if (loop == 1) { // discrete sol'n
                h += xt::linalg::dot(g, x0 - S.xc())();
            }
        } else { // can't cut in the previous iteration
            if (loop == 0) {
                break; // no more alternative cut
            }
            h += xt::linalg::dot(g, x0 - S.xc())();
        }
        if (t != t1) { // best t obtained
            feasible = true;
            t = t1;
            x_best = x0;
        }
        double tsq;
        std::tie(status, tsq) = S.update(g, h);
        if (status == 1) {
            break;
        }
        if (tsq < options.tol) {
            status = 2;
            break;
        }
    }

    return std::tuple{std::move(x_best), t, niter, feasible, status};
} // END

#endif