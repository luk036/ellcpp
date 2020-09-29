// -*- coding: utf-8 -*-
#pragma once

//#include "mat.hpp"
#include "chol_ext.hpp"
#include <boost/optional.hpp>
#include <gsl/span>
#include <xtensor/xarray.hpp>

/*!
 * @brief Oracle for Linear Matrix Inequality
 *
 *    This oracle solves the following feasibility problem:
 *
 *        find  x
 *        s.t.  F * x >= 0
 */
class lmi0_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Cut = std::tuple<Arr, double>;

  private:
    const gsl::span<const Arr> _F;
    const size_t _n;

  public:
    chol_ext<> _Q;

  public:
    /*!
     * @brief Construct a new lmi0 oracle object
     *
     * @param[in] F
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
     * @param[in] x
     * @return boost::optional<Cut>
     */
    auto operator()(const Arr& x) -> boost::optional<Cut>;
};
