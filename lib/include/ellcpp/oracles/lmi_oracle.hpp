#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI_ORACLE_HPP 1

//#include "mat.hpp"
#include "chol_ext.hpp"
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <vector>

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
    Arr &_F0;
    Arr _A;
    chol_ext _Q;

  public:
    explicit lmi_oracle(const std::vector<Arr> &F, Arr &B): 
      _F{F}, //
      _F0{B}, //
      _A{xt::zeros<double>(B.shape())}, //
      _Q(B.shape()[0]) // 
    {}

    auto operator()(const Arr &x) {
        using xt::linalg::dot;
        using xt::placeholders::_;
        auto n = x.size();

        auto getA = [&,this](unsigned i, unsigned j) -> double {
            this->_A(i, j) = this->_F0(i, j);
            for (auto k=0u; k < n; ++k) {
                const auto& Fi = _F[k];
                this->_A(i, j) -= Fi(i,j) * x(k);
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
            const auto& Fi = _F[i];
            g(i) = _Q.sym_quad(v, Fi);
        }
        return std::tuple{std::move(g), 1., false};
    }

    // auto chk_mtx(Arr A, const Arr &x) {
    //     using xt::linalg::dot;
    //     using xt::placeholders::_;

    //     auto n = x.size();
    //     Arr g = xt::zeros<double>({n});
    //     auto fj = -1.;
    //     for (auto i = 0u; i < n; ++i) {
    //         auto Fi = xt::view(_F, i, xt::all(), xt::all());
    //         // Arr Fi = _F(i);
    //         A -= Fi * x(i);
    //     }
    //     _Q.factorize(A);
    //     if (_Q.is_spd()) {
    //         return std::tuple{std::move(g), fj};
    //     }
    //     Arr v = _Q.witness();
    //     auto p = v.size();
    //     fj = 1.;
    //     for (auto i = 0u; i < n; ++i) {
    //         auto Fi = xt::view(_F, i, xt::all(), xt::all());
    //         g(i) = _Q.sym_quad(v, Fi);
    //     }
    //     return std::tuple{std::move(g), fj};
    // }

    // auto chk_spd_t(const Arr &x, double t) {
    //     Arr A = _F0;
    //     // ???
    //     // auto m = A.size();
    //     // for (auto i = 0u; i < m; ++i) {
    //     //   A(i, i) += t;
    //     // }
    //     A += t;
    //     return this->chk_mtx(A, x);
    // }

    // auto chk_spd(const Arr &x) { return this->chk_mtx(_F0, x); }

    // auto operator()(const Arr &x, double t) {
    //     auto [g, fj] = this->chk_spd_t(x, t);
    //     if (fj < 0) {
    //         t -= 1.;
    //     }
    //     return std::make_tuple(g, fj, t);
    // }
};

#endif