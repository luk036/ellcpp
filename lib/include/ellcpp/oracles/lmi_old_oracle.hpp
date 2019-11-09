// -*- coding: utf-8 -*-
#pragma once

#include "chol_ext.hpp"
#include <optional>
#include <xtensor/xarray.hpp>
#include <gsl/span>

/*!
 * @brief Oracle for Linear Matrix Inequality
 *
 * Oracle for:
 *    F * x <= B
 * or
 *    (B - F * x) must be a semidefinte matrix
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
    lmi_old_oracle(gsl::span<const Arr> F, Arr&& B)
        : _F {F}
        , _F0 {std::forward<Arr>(B)}
        , _Q {this->_F0.shape()[0]}
    {
    }

    /*!
     * @brief Construct a new lmi oracle object
     *
     * @param F
     * @param B
     */
    lmi_old_oracle(gsl::span<const Arr> F, const Arr& B)
        : _F {F}
        , _F0 {B}
        , _Q(this->_F0.shape()[0])
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
