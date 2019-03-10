#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI_ORACLE_HPP 1

#include "chol_ext.hpp"
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

/**
 * @brief  Oracle for Linear Matrix Inequality
 *
 * Oracle for:
 *    F * x <= B
 * or
 *    (B - F * x) must be a semidefinte matrix
 */
class lmi_oracle {
    using Arr = xt::xarray<double>;
    using shape_type = Arr::shape_type;

  private:
    const std::vector<Arr> &_F;
    Arr _F0;
    chol_ext _Q;

  public:
    /**
     * @brief Construct a new lmi oracle object
     *
     * @param F
     * @param B
     */
    lmi_oracle(const std::vector<Arr> &F, Arr &&B)
        : _F{F},             //
          _F0{std::move(B)}, //
          _Q(B.shape()[0])   //
    {}

    /**
     * @brief Construct a new lmi oracle object
     *
     * @param F
     * @param B
     */
    lmi_oracle(const std::vector<Arr> &F, const Arr &B)
        : _F{F},           //
          _F0{B},          //
          _Q(B.shape()[0]) //
    {}

    /**
     * @brief
     *
     * @param x
     * @return auto
     */
    auto operator()(const Arr &x) {
        using xt::linalg::dot;
        using xt::placeholders::_;
        auto n = x.size();

        auto getA = [&, this](unsigned i, unsigned j) -> double {
            auto a = this->_F0(i, j);
            for (auto k = 0U; k < n; ++k) {
                // const auto &Fi = _F[k];
                a -= this->_F[k](i, j) * x(k);
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
            // const auto &Fi = _F[i];
            g(i) = _Q.sym_quad(v, _F[i]);
        }
        return std::tuple{std::move(g), ep, false};
    }
};

#endif