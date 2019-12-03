// -*- coding: utf-8 -*-
#pragma once

#include <cmath>
#include <tuple>
#include <xtensor/xarray.hpp>

/*!
 * @brief Oracle for a profit maximization problem.
 *
 *    This example is taken from [Aliabadi and Salahi, 2013]:
 *
 *        max     p(A x1^α x2^β) − v1*x1 − v2*x2
 *        s.t.    x1 ≤ k
 *
 *    where:
 *
 *        p(A x1^α x2^β): Cobb-Douglas production function
 *        p: the market price per unit
 *        A: the scale of production
 *        α, β: the output elasticities
 *        x: input quantity
 *        v: output price
 *        k: a given constant that restricts the quantity of x1
 */
class profit_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Cut = std::tuple<Arr, double>;

  private:
    double _log_pA;
    double _log_k;
    Arr _v;

  public:
    Arr _a;

  public:
    /*!
     * @brief Construct a new profit oracle object
     *
     * @param p the market price per unit
     * @param A the scale of production
     * @param k a given constant that restricts the quantity of x1
     * @param a the output elasticities
     * @param v output price
     */
    profit_oracle(double p, double A, double k, Arr a, Arr v)
        : _log_pA {std::log(p * A)}
        , _log_k {std::log(k)}
        , _v {v}
        , _a {a}
    {
    }

    /*!
     * @brief
     *
     * @param y input quantity (in log scale)
     * @param t the best-so-far optimal value
     * @return std::tuple<Cut, double> Cut and the updated best-so-far value
     */
    auto operator()(const Arr& y, double t) const -> std::tuple<Cut, double>;
};

/*!
 * @brief Oracle for a profit maximization problem (robust version)
 *
 *    This example is taken from [Aliabadi and Salahi, 2013]:
 *
 *        max     p'(A x1^α' x2^β') - v1'*x1 - v2'*x2
 *        s.t.    x1 ≤ k'
 *
 *    where:
 *        α' = α ± e1
 *        β' = β ± e2
 *        p' = p ± e3
 *        k' = k ± e4
 *        v' = v ± e5
 *
 * @see profit_oracle
 */
class profit_rb_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

  private:
    Arr _uie;
    Arr _a;
    profit_oracle _P;

  public:
    /*!
     * @brief Construct a new profit rb oracle object
     *
     * @param p the market price per unit
     * @param A the scale of production
     * @param k a given constant that restricts the quantity of x1
     * @param a the output elasticities
     * @param v output price
     * @param ui paramters for uncertainty
     * @param e paramters for uncertainty
     * @param e3 paramters for uncertainty
     */
    profit_rb_oracle(double p, double A, double k, const Arr& a, const Arr& v,
        const Arr& e, double e3)
        : _uie {e}
        , _a {a}
        , _P(p - e3, A, k - e3, a, v + e3)
    {
    }

    /*!
     * @brief Make object callable for cutting_plane_dc()
     *
     * @param y input quantity (in log scale)
     * @param t the best-so-far optimal value
     * @return Cut and the updated best-so-far value
     *
     * @see cutting_plane_dc
     */
    auto operator()(const Arr& y, double t)
    {
        auto a_rb = this->_a;
        a_rb[0] += y[0] > 0. ? -this->_uie[0] : this->_uie[0];
        a_rb[1] += y[1] > 0. ? -this->_uie[1] : this->_uie[1];
        this->_P._a = a_rb;
        return this->_P(y, t);
    }
};

/*!
 * @brief Oracle for profit maximization problem (discrete version)
 *
 *    This example is taken from [Aliabadi and Salahi, 2013]
 *
 *        max     p(A x1^α x2^β) - v1*x1 - v2*x2
 *        s.t.    x1 ≤ k
 *
 *    where:
 *
 *        p(A x1^α x2^β): Cobb-Douglas production function
 *        p: the market price per unit
 *        A: the scale of production
 *        α, β: the output elasticities
 *        x: input quantity (must be integer value)
 *        v: output price
 *        k: a given constant that restricts the quantity of x1
 *
 * @see profit_oracle
 */
class profit_q_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Cut = std::tuple<Arr, double>;

  private:
    profit_oracle P;

  public:
    /*!
     * @brief Construct a new profit q oracle object
     *
     * @param p the market price per unit
     * @param A the scale of production
     * @param k a given constant that restricts the quantity of x1
     * @param a the output elasticities
     * @param v output price
     */
    profit_q_oracle(double p, double A, double k, const Arr& a, const Arr& v)
        : P(p, A, k, a, v)
    {
    }

    /*!
     * @brief Make object callable for cutting_plane_q()
     *
     * @param y input quantity (in log scale)
     * @param t the best-so-far optimal value
     * @return Cut and the updated best-so-far value
     *
     * @see cutting_plane_q
     */
    auto operator()(const Arr& y, double t, int /* unused */) const
        -> std::tuple<Cut, double, Arr, int>;
};
