#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI_ORACLE_HPP 1

//#include "mat.hpp"
#include <xtensor/xarray.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include "chol_ext.hpp"

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
  Arr& _F;
  Arr& _F0;

public:
  explicit lmi_oracle(Arr& F, Arr& B) : _F{F}, _F0{B} {}

  auto chk_mtx(Arr &A, const Arr &x) const {
    using xt::placeholders::_;
    using xt::linalg::dot;

    auto n = x.size();
    auto g = Arr(n);
    auto fj = -1.0;
    for (auto i = 0u; i < n; ++i) {
      auto Fi = xt::view(_F, xt::all(), xt::all(), i);
      A -= Fi * x[i];
    }
    chol_ext Q(A);
    if (Q.is_sd()) {
      return std::tuple{g, fj};
    }
    auto v = Q.witness();
    auto p = v.size();
    auto sub = xt::range(_, p);
    fj = -dot(v, dot(xt::view(A, sub, sub), v))();
    for (auto i = 0u; i < n; ++i) {
      auto Fi = xt::view(_F, xt::all(), xt::all(), i);
      g[i] = dot(v, dot(xt::view(Fi, sub, sub), v))();
    }
    return std::tuple{g, fj};
  }

  auto chk_spd_t(const Arr &x, double t) const {
    Arr A = _F0;
    // ???
    // auto m = A.size();
    // for (auto i = 0u; i < m; ++i) {
    //   A(i, i) += t;
    // }
    A += t;
    return this->chk_mtx(A, x);
  }

  auto chk_spd(const Arr &x) const {
    Arr A = _F0;
    return this->chk_mtx(A, x);
  }

  auto operator()(const Arr &x, double t) const {
    auto [g, fj] = this->chk_spd_t(x, t);
    if (fj < 0.0) {
      t -= 1.0;
    }
    return std::make_tuple(g, fj, t);
  }
};

#endif