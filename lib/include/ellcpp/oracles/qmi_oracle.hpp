// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_QMI_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_QMI_ORACLE_HPP 1

#include "chol_ext.hpp"
#include <cassert>
#include <vector>
#include <xtensor/xarray.hpp>

/**
 * @brief Oracle for Quadratic Matrix Inequality
 *
 *   F(x).T * F(x) <= I*t
 * where
 *   F(x) = F0 - (F1 * x1 + F2 * x2 + ...)
 */
class qmi_oracle {
    using Arr = xt::xarray<double>;
    using shape_type = Arr::shape_type;

  private:
    const std::vector<Arr> &_F;
    const Arr &_F0;
    double _t;
    std::size_t _nx = 0;
    std::size_t _count;
    Arr _Fx;
    // Arr _A;

  public:
    chol_ext _Q;

  public:
    /**
     * @brief Construct a new qmi oracle object
     *
     * @param F
     * @param F0
     */
    explicit qmi_oracle(const std::vector<Arr> &F, const Arr &F0)
        : _F{F}, _F0{F0}, _t{0.}, _count{0},
          _Q(F0.shape()[0]), _Fx{xt::zeros<double>(F0.shape())} {}

    /**
     * @brief Update t
     *
     * @param t
     */
    void update(double t) { _t = t; }

    /**
     * @brief
     *
     * @param x
     * @return auto
     */
    auto operator()(const Arr &x) -> std::tuple<Arr, double, bool>;
};

#endif
