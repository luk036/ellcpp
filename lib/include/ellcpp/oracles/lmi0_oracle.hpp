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
    const std::vector<Arr> &_F;
    size_t _n;

  public:
    Arr _A;
    chol_ext _Q;

  public:
    explicit lmi0_oracle(const std::vector<Arr> &F)
        : _F{F},                           //
          _n{F[0].shape()[0]},             //
          _A{xt::zeros<double>({_n, _n})}, //
          _Q(_n)                           //
    {}

    auto operator()(const Arr &x) {
        using xt::linalg::dot;
        using xt::placeholders::_;
        auto n = x.size();

        auto getA = [&, this](unsigned i, unsigned j) -> double {
            this->_A(i, j) = 0.;
            for (auto k = 0u; k < n; ++k) {
                // auto Fi = _F[k];
                this->_A(i, j) += _F[k](i, j) * x(k);
            }
            return this->_A(i, j);
        };

        Arr g = xt::zeros<double>({n});

        _Q.factor(getA);
        if (_Q.is_spd()) {
            return std::tuple{std::move(g), -1., true};
        }
        auto [v, ep] = _Q.witness();
        for (auto i = 0u; i < n; ++i) {
            // auto Fi = _F[i];
            g(i) = -_Q.sym_quad(v, _F[i]);
        }
        return std::tuple{std::move(g), ep, false};
    }
};

#endif