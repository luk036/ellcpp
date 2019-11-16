// -*- coding: utf-8 -*-
#pragma once

#include <cmath>
#include <tuple>
// #include <xtensor/xarray.hpp>

/*!
 * @brief Options
 *
 */
struct Options
{
    unsigned int max_it = 2000;
    double tol = 1e-8;
};

/*!
 * @brief CInfo
 *
 */
struct CInfo
{
    double value = 0.;
    bool feasible;
    size_t num_iters;
    int status;

    /*!
     * @brief Construct a new CInfo object
     *
     * @param feasible
     * @param num_iters
     * @param status
     */
    CInfo(bool feasible, size_t num_iters, int status)
        : feasible {feasible}
        , num_iters {num_iters}
        , status {status}
    {
    }
};

/*!
 * @brief
 *
 * @tparam Oracle
 * @tparam Space
 * @param Omega
 * @param I
 * @param options
 * @return auto
 */
template <typename Oracle, typename Space>
auto bsearch(Oracle& Omega, const Space& I, const Options& options = Options())
{
    // assume monotone
    auto [l, u] = I;
    const auto u_orig = u;
    auto niter = 1U;
    for (; niter <= options.max_it; ++niter)
    {
        auto t = l + (u - l) / 2;
        if (Omega(t))
        { // feasible sol'n obtained
            u = t;
        }
        else
        {
            l = t;
        }
        auto tau = (u - l) / 2;
        if (tau < options.tol)
        {
            break;
        }
    }
    auto ret = CInfo(u != u_orig, niter, 0);
    ret.value = u;
    return ret;
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
     * @param P
     * @param S
     * @param options
     */
    bsearch_adaptor(Oracle& P, Space& S, const Options& options = Options())
        : _P {P}
        , _S {S}
        , _options {options}
    {
    }

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
     * @return auto
     */
    auto operator()(double t)
    {
        Space S(this->_S);
        this->_P.update(t);
        auto ret_info = cutting_plane_feas(this->_P, S, this->_options);
        if (ret_info.feasible)
        {
            this->_S.set_xc(S.xc());
            return true;
        }
        return false;
    }
};

/*!
 * @brief Find a point in a convex set (defined through a cutting-plane oracle).
 *
 * @tparam Oracle
 * @tparam Space
 * @tparam T
 * @param[in] Omega   perform assessment on x0
 * @param[in] S        search Space containing x*
 * @param[in] max_it   maximum number of iterations
 * @param[in] tol      error tolerance
 * @return x      solution vector
 * @return niter  number of iterations performed
 * @return feasible   solution found or not
 * @return status how is the final cut
 */
template <typename Oracle, typename Space>
auto cutting_plane_feas(
    Oracle& Omega, Space& S, const Options& options = Options())
{
    auto feasible = false;
    auto niter = 1U, status = 0U;

    for (; niter <= options.max_it; ++niter)
    {
        auto cut = Omega(S.xc()); // query the oracle at S.xc()
        if (!cut)
        { // feasible sol'n obtained
            feasible = true;
            break;
        }
        double tsq;
        std::tie(status, tsq) = S.update(*cut); // update S
        if (status != 0)
        {
            break;
        }
        if (tsq < options.tol)
        { // no more
            status = 2;
            break;
        }
    }
    return CInfo(feasible, niter, status);
    // return std::tuple{S.xc(), niter, feasible, status};
}

/*!
 * @brief Cutting-plane method for solving convex problem
 *
 * @tparam Oracle
 * @tparam Space
 * @tparam T
 * @param[in] Omega    perform assessment on x0
 * @param[in] S        search Space containing x*
 * @param[in] t        best-so-far optimal sol'n
 * @param[in] max_it   maximum number of iterations
 * @param[in] tol      error tolerance
 * @return x      solution vector
 * @return niter   number of iterations performed
 * @return feasible   solution found or not
 * @return status how is the final cut
 */
template <typename Oracle, typename Space>
auto cutting_plane_dc(
    Oracle& Omega, Space& S, double t, const Options& options = Options())
{
    const auto t_orig = t;
    auto x_best = S.xc();
    auto niter = 1U;
    auto status = 0U;

    for (; niter <= options.max_it; ++niter)
    {
        auto [cut, t1] = Omega(S.xc(), t);
        if (t != t1)
        { // best t obtained
            // feasible = true;
            t = t1;
            x_best = S.xc();
        }
        double tsq;
        std::tie(status, tsq) = S.update(cut);
        if (status == 1)
        {
            break;
        }
        if (tsq < options.tol)
        { // no more
            status = 2;
            break;
        }
    }
    auto ret = CInfo(t != t_orig, niter, status);
    // ret.val = std::move(x_best);
    ret.value = t;
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
 * @brief
 *
 * @tparam Oracle
 * @tparam Space
 * @tparam T
 * @param Omega
 * @param S
 * @param t the best-so-far optimal value
 * @param options
 * @return auto
 */
template <typename Oracle, typename Space>
auto cutting_plane_q(
    Oracle& Omega, Space& S, double t, const Options& options = Options())
{
    const auto t_orig = t;
    auto x_best = S.xc();
    auto status = 1U;
    auto niter = 1U;

    for (; niter <= options.max_it; ++niter)
    {
        auto [cut, t1, x0, loop] = Omega(S.xc(), t, (status != 3) ? 0 : 1);
        // if (status != 3) {
        //     if (loop == 1) { // discrete sol'n
        //         h += xt::linalg::dot(g, x0 - S.xc())();
        //     }
        // } else { // can't cut in the previous iteration
        if (status == 3)
        {
            if (loop == 0)
            {
                break; // no more alternative cut
            }
            // h += xt::linalg::dot(g, x0 - S.xc())();
        }
        if (t != t1)
        { // best t obtained
            // feasible = true;
            t = t1;
            x_best = x0;
        }
        double tsq;
        std::tie(status, tsq) = S.update(cut);
        if (status == 1)
        {
            break;
        }
        if (tsq < options.tol)
        {
            status = 2;
            break;
        }
    }
    auto ret = CInfo(t != t_orig, niter, status);
    // ret.val = std::move(x_best);
    ret.value = t;
    return std::tuple {std::move(x_best), std::move(ret)};
} // END
