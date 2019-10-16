// -*- coding: utf-8 -*-
#pragma once

#include <cmath>
#include <tuple>
#include <xtensor/xarray.hpp>

/*!
 * @brief Ellipsoid Search Space
 *
 * ell = { x | (x - xc)' * (\kappa Q)^-1 * (x - xc) <= 1 }
 */
class ell
{
  public:
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using params_t = std::tuple<double, double, double>;
    using return_t = std::tuple<int, params_t>;

  public:
    bool _use_parallel_cut = true;
    bool _no_defer_trick = false;

  private:
    std::size_t _n;
    double _c1;
    double _kappa;
    Arr _xc;
    Arr _Q;

  public:
    /*!
     * @brief Construct a new ell object
     *
     * @tparam T
     * @param val
     * @param x
     */
    ell(const Arr& val, const Arr& x)
        : _n {x.size()}
        , _c1 {double(_n * _n) / (_n * _n - 1)}
        , _kappa {1}
        , _xc {x}
        , _Q {xt::diag(val)}
    {
    }

    /*!
     * @brief Construct a new ell object
     *
     * @tparam T
     * @param val
     * @param x
     */
    ell(const double& alpha, const Arr& x)
        : _n {x.size()}
        , _c1 {double(_n * _n) / (_n * _n - 1)}
        , _kappa {alpha}
        , _xc {x}
        , _Q {xt::eye(_n)}
    {
    }

    /*!
     * @brief Construct a new ell object
     *
     * @param E
     */
    ell(const ell& E) = default;

    /*!
     * @brief copy the whole array anyway
     *
     * @return const Arr&
     */
    auto xc() const -> Arr
    {
        return _xc;
    }

    // /*!
    //  * @brief
    //  *
    //  * @return Arr&
    //  */
    // auto xc() -> Arr&
    // {
    //     return _xc;
    // }

    /*!
     * @brief Set the xc object
     *
     * @param xc
     */
    void set_xc(const Arr& xc)
    {
        _xc = xc;
    }

    /*!
     * @brief Update ellipsoid core function using the cut
     *          g' * (x - xc) + beta <= 0
     *
     * @tparam T
     * @param g
     * @param beta
     * @return std::tuple<int, double>
     */
    template <typename T>
    std::tuple<int, double> update(const Arr& g, const T& beta);

  private:
    /*!
     * @brief
     *
     * @param b0
     * @param b1
     * @param tsq
     * @return return_t
     */
    return_t __calc_ll_core(double b0, double b1, double tsq) const;

    /*!
     * @brief Parallel Cut, one of them is central
     *
     * @param b1
     * @param t1
     * @param tsq
     * @return return_t
     */
    return_t __calc_ll_cc(double b1, double t1, double tsq) const;

    /*!
     * @brief Deep Cut
     *
     * @param b0
     * @param tsq
     * @return return_t
     */
    return_t __calc_dc(double b0, double tsq) const;

    /*!
     * @brief Central Cut
     *
     * @param tsq
     * @return return_t
     */
    return_t __calc_cc(double tsq) const;
}; // } ell
