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
    // using params_t = std::tuple<double, double, double>;
    // using return_t = std::tuple<int, params_t>;

  public:
    bool _use_parallel_cut = true;
    bool _no_defer_trick = false;

  private:
    double _rho = 0.;
    double _sigma = 0.;
    double _delta = 0.;
    double _tsq = 0.;

    size_t _n;
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
    ell(const Arr& val, Arr&& x)
        : _n {x.size()}
        , _c1 {double(_n * _n) / (_n * _n - 1)}
        , _kappa {1.}
        , _xc {std::forward<Arr>(x)}
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
    ell(const Arr& val, const Arr& x)
        : _n {x.size()}
        , _c1 {double(_n * _n) / (_n * _n - 1)}
        , _kappa {1.}
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
     * @tparam T
     * @param val
     * @param x
     */
    ell(const double& alpha, Arr&& x)
        : _n {x.size()}
        , _c1 {double(_n * _n) / (_n * _n - 1)}
        , _kappa {alpha}
        , _xc {std::forward<Arr>(x)}
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
    std::tuple<int, double> update(const std::tuple<Arr, T>& cut);

  private:
    /*!
     * @brief
     *
     * @param b0
     * @param b1
     * @return int
     */
    auto __calc_ll_core(const double& b0, const double& b1) -> int;

    /*!
     * @brief Parallel Cut, one of them is central
     *
     * @param b1
     * @param t1
     * @return void
     */
    auto __calc_ll_cc(const double& b1, const double& t1) -> void;

    /*!
     * @brief Deep Cut
     *
     * @param b0
     * @return int
     */
    auto __calc_dc(const double& b0) -> int;

    /*!
     * @brief Central Cut
     *
     * @param tau
     * @return void
     */
    auto __calc_cc(const double& tau) -> void;
}; // } ell
