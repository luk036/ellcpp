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
    chol_ext _Q;

  public:
    /**
     * @brief Construct a new lmi0 oracle object
     *
     * @param F
     */
    explicit lmi0_oracle(const std::vector<Arr> &F)
        : _F{F},               //
          _n{F[0].shape()[0]}, //
          _Q(_n)               //
    {}

    /**
     * @brief
     *
     * @param x
     * @return auto
     */
    auto operator()(const Arr &x) {
        auto n = x.size();

        auto getA = [&, this](unsigned i, unsigned j) -> double {
            auto a = 0.;
            for (auto k = 0U; k < n; ++k) {
                a += _F[k](i, j) * x(k);
            }
            return a;
        };

        Arr g = xt::zeros<double>({n});

        _Q.factor(getA);
        if (_Q.is_spd()) {
            return std::tuple{std::move(g), -1., true};
        }
        auto [v, ep] = _Q.witness();
        for (auto i = 0U; i < n; ++i) {
            // auto Fi = _F[i];
            g(i) = -_Q.sym_quad(v, _F[i]);
        }
        return std::tuple{std::move(g), ep, false};
    }
};

#endif