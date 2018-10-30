#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI0_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI0_ORACLE_HPP 1

//#include "mat.hpp"
#include "chol_ext.hpp"
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>


/**
 * @brief  Oracle for Linear Matrix Inequality
 *
 * Oracle for:
 *    F * x >= 0
 */
class lmi0_oracle {
    using Arr = xt::xarray<double>;
    using shape_type = Arr::shape_type;

  private:
    Arr &_F;
    size_t _n;
    Arr _A;
    chol_ext _Q;

  public:
    explicit lmi0_oracle(Arr &F) : 
      _F{F}, //
      _n{F.shape()[-1]}, //
      _A{xt::zeros<double>({_n, _n})}, //
      _Q(_n) // 
    {}

    auto operator()(const Arr &x) {
        using xt::linalg::dot;
        using xt::placeholders::_;
        auto n = x.size();

        auto getA = [&,this](unsigned i, unsigned j) -> double {
            this->_A(i, j) = 0.;
            for (auto k=0u; k < n; ++k) {
                auto Fi = xt::view(_F, k, xt::all(), xt::all());
                this->_A(i, j) += Fi(i,j) * x(k);
            }
            return this->_A(i, j);
        };

        Arr g = xt::zeros<double>({n});

        _Q.factor(getA);
        if (_Q.is_spd()) {
            return std::tuple{std::move(g), -1., true};
        }
        Arr v = _Q.witness();
        for (auto i = 0u; i < n; ++i) {
            auto Fi = xt::view(_F, i, xt::all(), xt::all());
            g(i) = -_Q.sym_quad(v, Fi);
        }
        return std::tuple{std::move(g), 1., false};
    }
};

#endif