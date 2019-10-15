// -*- coding: utf-8 -*-
#pragma once

#include "chol_ext.hpp"
#include <vector>
#include <xtensor/xarray.hpp>

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

  private:
    const std::vector<Arr>& _F;
    Arr _F0;
    chol_ext<> _Q;

  public:
    /*!
     * @brief Construct a new lmi oracle object
     *
     * @param F
     * @param B
     */
    lmi_old_oracle(const std::vector<Arr>& F, Arr&& B)
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
    lmi_old_oracle(const std::vector<Arr>& F, const Arr& B)
        : _F {F}
        , _F0 {B}
        , _Q(this->_F0.shape()[0])
    {
    }

    /*!
     * @brief
     *
     * @param x
     * @return auto
     */
    auto operator()(const Arr& x) -> std::tuple<bool, Arr, double>;
};
