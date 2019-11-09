// -*- coding: utf-8 -*-
#pragma once

//#include "mat.hpp"
#include "chol_ext.hpp"
#include <optional>
#include <xtensor/xarray.hpp>
#include <gsl/span>

/*!
 * @brief  Oracle for Linear Matrix Inequality
 *
 *   find  x
 *   s.t.​  F*x ⪰ 0
 */
class lmi0_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Cut = std::tuple<Arr, double>;

  private:
    gsl::span<const Arr> _F;
    size_t _n;

  public:
    chol_ext<> _Q;

  public:
    /*!
     * @brief Construct a new lmi0 oracle object
     *
     * @param F
     */
    explicit lmi0_oracle(gsl::span<const Arr> F)
        : _F {F}
        , _n {F[0].shape()[0]}
        , _Q(_n)
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
