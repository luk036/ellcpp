// -*- coding: utf-8 -*-
#pragma once

#include "chol_ext.hpp"
#include <gsl/span>
#include <optional>
#include <xtensor/xarray.hpp>

/*!
 * @brief Oracle for Linear Matrix Inequality.
 *
 *    This oracle solves the following feasibility problem:
 *
 *        find  x
 *        s.t.  (B − F * x) ⪰ 0
 */
class lmi_old_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Cut = std::tuple<Arr, double>;

  private:
    gsl::span<const Arr> _F;
    Arr _F0;
    chol_ext<> _Q;

  public:
    /*!
     * @brief Construct a new lmi oracle object
     *
     * @param F
     * @param B
     */
    lmi_old_oracle(gsl::span<const Arr> F, Arr B)
        : _F {F}
        , _F0 {std::move(B)}
        , _Q {this->_F0.shape()[0]}
    {
    }

    /*!
     * @brief
     *
     * @param x
     * @return std::optional<Cut>
     */
    auto operator()(const Arr& x) -> std::optional<Cut>;
};
