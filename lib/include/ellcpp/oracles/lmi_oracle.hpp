// -*- coding: utf-8 -*-
#pragma once

#include "chol_ext.hpp"
#include <gsl/span>
#include <boost/optional.hpp>
#include <xtensor/xarray.hpp>

/*!
 * @brief Oracle for Linear Matrix Inequality.
 *
 *    This oracle solves the following feasibility problem:
 *
 *        find  x
 *        s.t.  (B - F * x) >= 0
 */
class lmi_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Cut = std::tuple<Arr, double>;

  private:
    const gsl::span<const Arr> _F;
    const Arr _F0;
    chol_ext<> _Q;

  public:
    /*!
     * @brief Construct a new lmi oracle object
     *
     * @param[in] F
     * @param[in] B
     */
    lmi_oracle(gsl::span<const Arr> F, Arr B)
        : _F {F}
        , _F0 {std::move(B)}
        , _Q {this->_F0.shape()[0]}
    {
    }

    /*!
     * @brief
     *
     * @param[in] x
     * @return boost::optional<Cut>
     */
    auto operator()(const Arr& x) -> boost::optional<Cut>;
};
