// -*- coding: utf-8 -*-
#pragma once

#include <cmath>
#include <tuple>
#include <xtensor/xarray.hpp>

// forward declaration
enum class CUTStatus;

/*!
 * @brief Ellipsoid Search Space
 *
 *        ell = {x | (x - xc)' Q^-1 (x - xc) \le \kappa}
 */
class ell
{
  public:
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    // using params_t = std::tuple<double, double, double>;
    // using return_t = std::tuple<int, params_t>;

    bool use_parallel_cut = true;
    bool no_defer_trick = false;

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
     * @param[in] val
     * @param[in] x
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
     * @param[in] alpha
     * @param[in] x
     */
    ell(const double& alpha, Arr x)
        : _n {x.size()}
        , _c1 {double(_n * _n) / (_n * _n - 1)}
        , _kappa {alpha}
        , _xc {std::move(x)}
        , _Q {xt::eye(_n)}
    {
    }

    /**
     * @brief Construct a new ell object
     *
     * @param[in] E (move)
     */
    ell(ell&& E) = default;

  private:
    /*!
     * @brief Construct a new ell object
     *
     * @param[in] E
     */
    ell(const ell& E) = default;
    ell& operator=(const ell& E) = delete;

  public:
    /**
     * @brief explicitly copy
     *
     * @return ell
     */
    ell copy() const
    {
        return ell(*this);
    }

    /*!
     * @brief copy the whole array anyway
     *
     * @return Arr
     */
    auto xc() const -> Arr
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
    std::tuple<CUTStatus, double> update(const std::tuple<Arr, T>& cut);

  private:
    /*!
     * @brief Calculate new ellipsoid under Parallel Cut
     *
     *        g' (x - xc​) + beta0 \le 0
     *        g' (x - xc​) + beta1 \ge 0
     *
     * @param[in] b0
     * @param[in] b1
     * @return int
     */
    auto _calc_ll_core(const double& b0, const double& b1) -> CUTStatus;

    /*!
     * @brief Calculate new ellipsoid under Parallel Cut, one of them is central
     *
     *        g' (x - xc​) \le 0
     *        g' (x - xc​) + beta1 \ge 0
     *
     * @param[in] b1
     * @param[in] b1sq
     */
    auto _calc_ll_cc(const double& b1, const double& b1sq) -> void;

    /*!
     * @brief Calculate new ellipsoid under Deep Cut
     *
     *        g' (x - xc​) + beta \le 0
     *
     * @param[in] beta
     */
    auto _calc_dc(const double& beta) -> CUTStatus;

    /*!
     * @brief Calculate new ellipsoid under Central Cut
     *
     *        g' (x - xc​) \le 0
     *
     * @param[in] tau
     */
    auto _calc_cc(const double& tau) -> void;
}; // } ell
