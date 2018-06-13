#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_CUTTING_PLANE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_CUTTING_PLANE_HPP 1

#include <cmath>
#include <tuple>

struct Options {
    unsigned int max_it = 2000;
    double tol = 1e-4;
};

template <typename Oracle, typename Space>
auto bsearch(Oracle &evaluate, const Space &I,
             const Options &options = Options()) {
    // assume monotone
    bool feasible = false;
    auto [l, u] = I;
    auto t = l + (u - l) / 2;
    auto niter = 1u;
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

template <typename Oracle, typename Space> class bsearch_adaptor {
  private:
    Oracle &_P;
    Space &_S;
    Options _options;

  public:
    explicit bsearch_adaptor(Oracle &P, Space &S,
                             const Options &options = Options())
        : _P{P}, _S{S}, _options{options} {}

    auto x_best() const { return _S.xc(); }

    auto operator()(double t) {
        Space S(_S);
        _P.update(t);
        auto [x, _1, feasible, _2] = cutting_plane_feas(_P, S, _options);
        if (feasible) {
            _S.set_xc(x);
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
    auto niter = 1u, status = 0u;

    for (; niter <= options.max_it; ++niter) {
        auto [g, h, flag] = evaluate(S.xc());
        if (flag) { // feasible sol'n obtained
            feasible = true;
            break;
        }
        double tau;
        std::tie(status, tau) = S.update(g, h);
        if (status != 0) {
            break;
        }
        if (tau < options.tol) { // no more
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
    bool feasible = false;
    auto x_best = S.xc();
    auto niter = 1u, status = 0u;
    for (; niter <= options.max_it; ++niter) {
        auto [g, h, t1] = evaluate(S.xc(), t);
        if (t != t1) { // best t obtained
            feasible = true;
            t = t1;
            x_best = S.xc();
        }
        double tau;
        std::tie(status, tau) = S.update(g, h);
        if (status == 1) {
            break;
        }
        if (tau < options.tol) { // no more
            status = 2;
            break;
        }
    }
    return std::tuple{x_best, t, niter, feasible, status};
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
#include <xtensor/xarray.hpp>

template <typename Oracle, typename Space, typename T>
auto cutting_plane_q(Oracle &evaluate, Space &S, T t,
                     const Options &options = Options()) {
    bool feasible = false;
    auto x_best = S.xc();
    auto niter = 1u, status = 1u;
    for (; niter < options.max_it; ++niter) {
        auto [g, h, t1, x, loop] = evaluate(S.xc(), t, (status != 3) ? 0 : 1);
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
            feasible = true;
            t = t1;
            x_best = x;
        }
        double tau;
        std::tie(status, tau) = S.update(g, h);
        if (status == 1) {
            break;
        }
        if (tau < options.tol) {
            status = 2;
            break;
        }
    }

    return std::tuple{x_best, t, niter, feasible, status};
} // END

#endif