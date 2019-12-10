// -*- coding: utf-8 -*-
#pragma once

#include <cmath>
#include <tuple>
#include <xtensor/xarray.hpp>

/*!
 * @brief Ellipsoid Search Space
 *
 *        ell = {x | (x − xc)' Q^−1 (x − xc) ​≤ κ}
 */
class ell
{
  public:
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    // using params_t = std::tuple<double, double, double>;
    // using return_t = std::tuple<int, params_t>;

    bool _use_parallel_cut = true;
    bool _no_defer_trick = false;

  private:
    double _rho {};
    double _sigma {};
    double _delta {};
    double _tsq {};

    const size_t _n;
    const double _c1;
    double _kappa;
    Arr _xc;
    Arr _Q;

  public:
    /*!
     * @brief Construct a new ell object
     *
     * @param val
     * @param x
     */
    ell(const Arr& val, Arr x)
        : _n {x.size()}
        , _c1 {double(_n * _n) / (_n * _n - 1)}
        , _kappa {1.}
        , _xc {std::move(x)}
        , _Q {xt::diag(val)}
    {
    }

    /*!
     * @brief Construct a new ell object
     *
     * @param alpha
     * @param x
     */
    ell(const double& alpha, Arr x)
        : _n {x.size()}
        , _c1 {double(_n * _n) / (_n * _n - 1)}
        , _kappa {alpha}
        , _xc {std::move(x)}
        , _Q {xt::eye(_n)}
    {
    }


  private:
    /*!
     * @brief Construct a new ell object
     *
     * @param E
     */
    ell(const ell& E) = default;
    ell& operator=(const ell& E) = delete;

  public:
    ell copy() const { return ell(*this); }

    /*!
     * @brief copy the whole array anyway
     *
     * @return const Arr&
     */
    auto xc() const -> Arr
    {
        return _xc;
    }

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
     * @brief Update ellipsoid core function using the cut(s)
     *
     * @tparam T
     * @param cut cutting-plane
     * @return std::tuple<int, double>
     */
    template <typename T>
    std::tuple<int, double> update(const std::tuple<Arr, T>& cut);

  private:
    /*!
     * @brief Calculate new ellipsoid under Parallel Cut
     *
     *        g' (x − xc​) + β0 ​≤ 0
     *        g' (x − xc​) + β1 ​≥ 0
     *
     * @param b0
     * @param b1
     * @return int
     */
    auto __calc_ll_core(const double& b0, const double& b1) -> int;

    /*!
     * @brief Calculate new ellipsoid under Parallel Cut, one of them is central
     *
     *        g' (x − xc​) ​≤ 0
     *        g' (x − xc​) + β1 ​≥ 0
     *
     * @param b1
     * @param b1sq
     */
    auto __calc_ll_cc(const double& b1, const double& b1sq) -> void;

    /*!
     * @brief Calculate new ellipsoid under Deep Cut
     *
     *        g' (x − xc​) + β ​≤ 0
     *
     * @param beta
     */
    auto __calc_dc(const double& beta) -> int;

    /*!
     * @brief Calculate new ellipsoid under Central Cut
     *
     *        g' (x − xc​) ≤ 0
     *
     * @param tau
     */
    auto __calc_cc(const double& tau) -> void;
}; // } ell
