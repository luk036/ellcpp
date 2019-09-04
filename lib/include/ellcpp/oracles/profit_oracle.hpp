// -*- coding: utf-8 -*-
#pragma once

//#include <valarray>
#include <cmath>
#include <tuple>
#include <xtensor/xarray.hpp>

/*!
 * @brief Oracle for a profit maximization problem
 *
 */
class profit_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

private:
    double _log_pA;
    double _log_k;
    Arr    _v;

public:
    Arr _a;

public:
    /*!
     * @brief Construct a new profit oracle object
     *
     * @param p
     * @param A
     * @param k
     * @param a
     * @param v
     */
    profit_oracle(double p, double A, double k, Arr a, Arr v)
        : _log_pA{std::log(p * A)}, //
          _log_k{std::log(k)},      //
          _v{v},                    //
          _a{a}                     //
    {
    }

    /*!
     * @brief
     *
     * @param y
     * @param t
     * @return auto
     */
    auto operator()(const Arr& y, double t) const -> std::tuple<Arr, double, double>;
};

/*!
 * @brief Oracle for a profit maximization problem (robust version)
 *
 */
class profit_rb_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

private:
    Arr    _uie;
    Arr    _a;
    // double _uie3;

    profit_oracle _P;

public:
    /*!
     * @brief Construct a new profit rb oracle object
     *
     * @param p
     * @param A
     * @param k
     * @param a
     * @param v
     * @param ui
     * @param e
     * @param e3
     */
    profit_rb_oracle(double p, double A, double k, const Arr& a, const Arr& v, const Arr& e,
                     double e3)
        : _uie{e},                         //
          _a{a},                           //
          _P(p - e3, A, k - e3, a, v + e3) //
    {
    }

    /*!
     * @brief
     *
     * @param y
     * @param t
     * @return auto
     */
    auto operator()(const Arr& y, double t)
    {
        auto a_rb = _a;
        for (auto&& i : {0, 1})
        {
            a_rb[i] += y[i] > 0 ? -_uie[i] : _uie[i];
        }
        _P._a = a_rb;
        return _P(y, t);
    }
};

/*!
 * @brief Oracle for profit maximization problem (discrete version)
 *
 */
class profit_q_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

private:
    profit_oracle P;

public:
    /*!
     * @brief Construct a new profit q oracle object
     *
     * @param p
     * @param A
     * @param k
     * @param a
     * @param v
     */
    profit_q_oracle(double p, double A, double k, const Arr& a, const Arr& v) : P(p, A, k, a, v) {}

    /*!
     * @brief
     *
     * @param y
     * @param t
     * @return auto
     */
    auto operator()(const Arr& y, double t, int /*unused*/) const
        -> std::tuple<Arr, double, double, Arr, int>;
};
