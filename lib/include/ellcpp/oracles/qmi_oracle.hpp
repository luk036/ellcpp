// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_QMI_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_QMI_ORACLE_HPP 1

#include "chol_ext.hpp"
#include <cassert>
// #include <vector>
#include <xtensor-blas/xlinalg.hpp>
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
    const Arr &_F;
    const Arr &_F0;
    double _t;
    std::size_t _nx = 0;
    std::size_t _count;
    Arr _Fx;
    Arr _A;

  public:
    chol_ext _Q;

  public:
    explicit qmi_oracle(const Arr &F, const Arr &F0)
        : _F{F}, _F0{F0}, _t{0.}, _count{0}, _A{xt::zeros<double>(F0.shape())},
          _Q(F0.shape()[0]), _Fx{xt::zeros<double>(F0.shape())} {}

    void update(double t) { _t = t; }

    auto operator()(const Arr &x) {
        using xt::linalg::dot;
        using xt::placeholders::_;

        _count = 0;
        _nx = x.shape()[0];

        auto getA = [this, x](std::size_t i, std::size_t j) -> double { // ???
            using xt::linalg::dot;
            assert(i >= j);
            Arr Fxi = xt::view(_Fx, i, xt::all());
            Arr Fxj = xt::view(_Fx, j, xt::all());
            if (_count < i + 1) {
                _count = i + 1;
                Fxi = xt::view(_F0, i, xt::all());
                for (auto k = 0u; k < _nx; ++k) {
                    // Arr Fk = _F[k];
                    Arr Fki = xt::view(_F, k, i, xt::all());
                    Fxi -= Fki * x(k);
                }
            }
            _A(i, j) = -dot(Fxi, Fxj)();
            if (i == j) {
                _A(i, j) += _t;
            }
            return _A(i, j);
        };

        _Q.factor(getA);
        Arr g = xt::zeros<double>({_nx});

        if (_Q.is_spd()) {
            return std::tuple{g, -1., true};
        }
        Arr v = _Q.witness();
        auto p = v.size();
        Arr Fxp = xt::transpose(xt::view(_Fx, xt::range(_, p), xt::all()));
        Arr Av = dot(Fxp, v);
        for (auto k = 0u; k < _nx; ++k) {
            // Arr Fk = _F[k];
            Arr Fkp = xt::view(_F, k, xt::range(_, p), xt::all());
            g(k) = -2. * dot(v, dot(Fkp, Av))();
        }
        return std::tuple{g, 1., false};
    }
};

#endif
