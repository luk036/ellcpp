// -*- coding: utf-8 -*-
#pragma once

#include <cmath>
#include <ellcpp/utility.hpp>
#include <tuple>
#include <xtensor/xarray.hpp>

// forward declaration
enum class CUTStatus;

/*!
 * @brief Ellipsoid Search Space
 *
 *        ell = {x | (x - xc)' Q^-1 (x - xc) \le \kappa}
 *
 * Keep $Q$ symmetric but no promise of positive definite
 */
class ell
{
  public:
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    // using params_t = std::tuple<double, double, double>;
    // using return_t = std::tuple<int, params_t>;

    bool use_parallel_cut = true;
    bool no_defer_trick = false;

  protected:
    double _mu {};
    double _rho {};
    double _sigma {};
    double _delta {};
    double _tsq {};

    const int _n;

    const double _nFloat;
    const double _nPlus1;
    const double _nMinus1;
    const double _halfN;
    const double _halfNplus1;
    const double _halfNminus1;
    const double _nSq;
    const double _c1;
    const double _c2;
    const double _c3;

    double _kappa;
    Arr _Q;
    Arr _xc;

    /*!
     * @brief Construct a new ell object
     *
     * @param[in] E
     */
    auto operator=(const ell& E) -> ell& = delete;

    /*!
     * @brief Construct a new ell object
     *
     * @param[in] val
     * @param[in] x
     */
    template <typename V, typename U>
    ell(V&& kappa, Arr&& Q, U&& x) noexcept
        : _n {int(x.size())}
        , _nFloat {double(_n)}
        , _nPlus1 {_nFloat + 1.}
        , _nMinus1 {_nFloat - 1.}
        , _halfN {_nFloat / 2.}
        , _halfNplus1 {_nPlus1 / 2.}
        , _halfNminus1 {_nMinus1 / 2.}
        , _nSq {_nFloat * _nFloat}
        , _c1 {_nSq / (_nSq - 1)}
        , _c2 {2. / _nPlus1}
        , _c3 {_nFloat / _nPlus1}
        , _kappa {std::forward<V>(kappa)}
        , _Q {std::move(Q)}
        , _xc {std::forward<U>(x)}
    {
    }

  public:
    /*!
     * @brief Construct a new ell object
     *
     * @param[in] val
     * @param[in] x
     */
    ell(const Arr& val, Arr x) noexcept
        : ell {1., xt::diag(val), std::move(x)}
    {
    }

    /*!
     * @brief Construct a new ell object
     *
     * @param[in] alpha
     * @param[in] x
     */
    ell(const double& alpha, Arr x) noexcept
        : ell {alpha, xt::eye(x.size()), std::move(x)}
    {
    }

    /**
     * @brief Construct a new ell object
     *
     * @param[in] E (move)
     */
    ell(ell&& E) = default;

    /**
     * @brief Destroy the ell object
     *
     */
    ~ell() { }

    /**
     * @brief Construct a new ell object
     *
     * To avoid accidentally copying, only explicit copy is allowed
     *
     * @param E
     */
    explicit ell(const ell& E) = default;

    /**
     * @brief explicitly copy
     *
     * @return ell
     */
    [[nodiscard]] auto copy() const -> ell
    {
        return ell(*this);
    }

    /*!
     * @brief copy the whole array anyway
     *
     * @return Arr
     */
    [[nodiscard]] auto xc() const -> Arr
    {
        return _xc;
    }

    /*!
     * @brief Set the xc object
     *
     * @param[in] xc
     */
    void set_xc(const Arr& xc)
    {
        _xc = xc;
    }

    /*!
     * @brief Update ellipsoid core function using the cut(s)
     *
     * @tparam T
     * @param[in] cut cutting-plane
     * @return std::tuple<int, double>
     */
    template <typename T>
    auto update(const std::tuple<Arr, T>& cut) -> std::tuple<CUTStatus, double>;

  protected:
    auto _update_cut(const double& beta) -> CUTStatus
    {
        return this->_calc_dc(beta);
    }

    CUTStatus _update_cut(const Arr& beta)
    { // parallel cut
        if (beta.shape()[0] < 2)
        {
            return this->_calc_dc(beta[0]);
        }
        return this->_calc_ll_core(beta[0], beta[1]);
    }

    /*!
     * @brief Calculate new ellipsoid under Parallel Cut
     *
     *        g' (x - xc) + beta0 \le 0
     *        g' (x - xc) + beta1 \ge 0
     *
     * @param[in] b0
     * @param[in] b1
     * @return int
     */
    auto _calc_ll_core(const double& b0, const double& b1) -> CUTStatus;

    /*!
     * @brief Calculate new ellipsoid under Parallel Cut, one of them is central
     *
     *        g' (x - xc) \le 0
     *        g' (x - xc) + beta1 \ge 0
     *
     * @param[in] b1
     * @param[in] b1sq
     */
    auto _calc_ll_cc(const double& b1, const double& b1sq) -> void;

    /*!
     * @brief Calculate new ellipsoid under Deep Cut
     *
     *        g' (x - xc) + beta \le 0
     *
     * @param[in] beta
     */
    auto _calc_dc(const double& beta) noexcept -> CUTStatus;

    /*!
     * @brief Calculate new ellipsoid under Central Cut
     *
     *        g' (x - xc) \le 0
     *
     * @param[in] tau
     */
    auto _calc_cc(const double& tau) noexcept -> void;
}; // } ell
