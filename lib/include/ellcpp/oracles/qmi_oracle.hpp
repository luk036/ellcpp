// -*- coding: utf-8 -*-
#pragma once

#include "ldlt_ext.hpp"
#include <gsl/span>
#include <optional>
#include <xtensor/xarray.hpp>

/*!
 * @brief Oracle for Quadratic Matrix Inequality
 *
 *    This oracle solves the following feasibility problem:
 *
 *        find  x
 *        s.t.  t * I - F(x)' F(x) >= 0
 *
 *    where
 *
 *        F(x) = F0 - (F1 * x1 + F2 * x2 + ...)
 */
class qmi_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Cut = std::tuple<Arr, double>;

  private:
    double _t = 0.;
    size_t _nx = 0;
    size_t _count = 0;

    const size_t _n;
    const size_t _m;
    const gsl::span<const Arr> _F;
    const Arr _F0;
    Arr _Fx;

  public:
    ldlt_ext _Q;


    /*!
     * @brief Construct a new qmi oracle object
     *
     * @param[in] F
     * @param[in] F0
     */
    qmi_oracle(gsl::span<const Arr> F, Arr F0)
        : _n {F0.shape()[0]}
        , _m {F0.shape()[1]}
        , _F {F}
        , _F0 {std::move(F0)}
        , _Fx {zeros({_m, _n})} // transposed
        , _Q(_m)                // take column
    {
    }

    /*!
     * @brief Update t
     *
     * @param[in] t the best-so-far optimal value
     */
    auto update(double t) -> void
    {
        this->_t = t;
    }

    /*!
     * @brief
     *
     * @param[in] x
     * @return std::optional<Cut>
     */
    auto operator()(const Arr& x) -> std::optional<Cut>;
};
