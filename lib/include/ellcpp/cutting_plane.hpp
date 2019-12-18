// -*- coding: utf-8 -*-
#pragma once

#include <cassert>
#include <cmath>
#include <tuple>
#include "half_nonnegative.hpp"

/*!
 * @brief Options
 *
 */
struct Options
{
    unsigned int max_it = 2000; //!< maximum number of iterations
    double tol = 1e-8;          //!< error tolerance
};

/*!
 * @brief CInfo
 *
 */
struct CInfo
{
    bool feasible;
    size_t num_iters;
    unsigned int status;
};

/*!
 * @brief
 *
 * @tparam Oracle
 * @tparam Space
 * @param Omega
 * @param I
 * @param options
 * @return CInfo
 */
template <typename Oracle, typename Space>
auto bsearch(Oracle&& Omega, Space&& I, const Options& options = Options())
    -> CInfo
{
    // assume monotone
    auto& [lower, upper] = I;
    assert(lower <= upper);
    const auto u_orig = upper;
    auto niter = 0U;
    auto status = 0U;

    for (; niter != options.max_it; ++niter)
    {
        auto tau = algo::half_nonnegative(upper - lower);
        if (tau < options.tol)
        {
            status = 2;
            break;
        }

        auto t = lower; // l may be `int` or `Fraction`
        t += tau;
        if (Omega(t))
        { // feasible sol'n obtained
            upper = t;
        }
        else
        {
            lower = t;
        }
    }
    return {upper != u_orig, niter + 1, status};
}

/*!
 * @brief
 *
 * @tparam Oracle
 * @tparam Space
 */
template <typename Oracle, typename Space> //
class bsearch_adaptor
{
  private:
    Oracle& _P;
    Space& _S;
    Options _options;

  public:
    /*!
     * @brief Construct a new bsearch adaptor object
     *
     * @param P perform assessment on x0
     * @param S search Space containing x*
     * @param options
     */
    bsearch_adaptor(Oracle& P, Space& S, const Options& options = Options())
        : _P {P}
        , _S {S}
        , _options {options}
    {
    }

    bsearch_adaptor(const bsearch_adaptor<Oracle, Space>& ) = delete;
    bsearch_adaptor& operator=(const bsearch_adaptor<Oracle, Space>& ) = delete;

    /*!
     * @brief get best x
     *
     * @return auto
     */
    auto x_best() const
    {
        return this->_S.xc();
    }

    /*!
     * @brief
     *
     * @param t the best-so-far optimal value
     * @return bool
     */
    template <typename opt_type>
    auto operator()(const opt_type& t) -> bool
    {
        Space S = this->_S.copy();
        this->_P.update(t);
        const auto ell_info = cutting_plane_feas(this->_P, S, this->_options);
        if (ell_info.feasible)
        {
            this->_S.set_xc(S.xc());
        }
        return ell_info.feasible;
    }
};

/*!
 * @brief Find a point in a convex set (defined through a cutting-plane oracle).
 * 
 *     A function f(x) is *convex* if there always exist a g(x)
 *     such that f(z) >= f(x) + g(x)' * (z - x), forall z, x in dom f.
 *     Note that dom f does not need to be a convex set in our definition.
 *     The affine function g' (x - xc) + beta is called a cutting-plane,
 *     or a ``cut'' for short.
 *     This algorithm solves the following feasibility problem:
 *
 *             find x
 *             s.t. f(x) <= 0,
 *
 *     A *separation oracle* asserts that an evalution point x0 is feasible,
 *     or provide a cut that separates the feasible region and x0.
 *    
 * @tparam Oracle 
 * @tparam Space 
 * @param Omega    perform assessment on x0
 * @param S        search Space containing x*
 * @param options  Maximum iteration and error tolerance etc.
 * @return Information of Cutting-plane method
 */
template <typename Oracle, typename Space>
auto cutting_plane_feas(
    Oracle&& Omega, Space&& S, const Options& options = Options()) -> CInfo
{
    auto feasible = false;
    auto niter = 1U;
    auto status = 0U;

    for (; niter != options.max_it; ++niter)
    {
        auto cut = Omega(S.xc()); // query the oracle at S.xc()
        if (not cut)
        { // feasible sol'n obtained
            feasible = true;
            break;
        }
        const auto [cutstatus, tsq] = S.update(*cut); // update S
        if (cutstatus != 0)
        {
            status = cutstatus;
            break;
        }
        if (tsq < options.tol)
        { // no more
            status = 2;
            break;
        }
    }
    return {feasible, niter, status};
}

/*!
 * @brief Cutting-plane method for solving convex problem
 *
 * @tparam Oracle
 * @tparam Space
 * @tparam opt_type
 * @param Omega    perform assessment on x0
 * @param S        search Space containing x*
 * @param t[inout] best-so-far optimal sol'n
 * @param options  Maximum iteration and error tolerance etc.
 * @return Information of Cutting-plane method
 */
template <typename Oracle, typename Space, typename opt_type>
auto cutting_plane_dc(
    Oracle&& Omega, Space&& S, opt_type&& t, const Options& options = Options())
{
    const auto t_orig = t;
    auto x_best = S.xc();
    auto niter = 1U;
    auto status = 0U;

    for (; niter != options.max_it; ++niter)
    {
        const auto [cut, t1] = Omega(S.xc(), t);
        if (t != t1)
        { // best t obtained
            // feasible = true;
            t = t1;
            x_best = S.xc();
        }
        const auto [cutstatus, tsq] = S.update(cut);
        if (cutstatus != 0) // ???
        {
            status = cutstatus;
            break;
        }
        if (tsq < options.tol)
        { // no more
            status = 2;
            break;
        }
    }
    auto ret = CInfo{t != t_orig, niter, status};
    return std::tuple {std::move(x_best), std::move(ret)};
} // END

/*!
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
// #include <xtensor-blas/xlinalg.hpp>
// #include <xtensor/xarray.hpp>

/*!
 * @brief Cutting-plane method for solving convex discrete optimization problem
 *
 * @tparam Oracle
 * @tparam Space
 * @param Omega    perform assessment on x0
 * @param S        search Space containing x*
 * @param t[inout] best-so-far optimal sol'n
 * @param options  Maximum iteration and error tolerance etc.
 * @return Information of Cutting-plane method
 */
template <typename Oracle, typename Space, typename opt_type>
auto cutting_plane_q(
    Oracle&& Omega, Space&& S, opt_type&& t, const Options& options = Options())
{
    auto x_best = S.xc();  // copying
    const auto t_orig = t;
    auto status = 1U; // note!!!
    auto niter = 1U;

    for (; niter != options.max_it; ++niter)
    {
        const auto [cut, x0, t1, loop] =
            Omega(S.xc(), t, (status != 3) ? 0 : 1);
        if (status == 3)
        {
            if (loop == 0)
            {
                break; // no more alternative cut
            }
        }
        if (t != t1)
        { // best t obtained
            t = t1;
            x_best = x0;
        }
        const auto [cutstatus, tsq] = S.update(cut);
        if (cutstatus == 1)
        {
            status = cutstatus;
            break;
        }
        if (tsq < options.tol)
        {
            status = 2;
            break;
        }
    }
    auto ret = CInfo{t != t_orig, niter, status};
    return std::tuple {std::move(x_best), std::move(ret)};
} // END
