#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI_ORACLE_HPP 1

//#include "mat.hpp"
#include "chol_ext.hpp"
#include <xtensor/xarray.hpp>
#include <xtensor-blas/xlinalg.hpp>

class lmi_oracle {
  using Vec = xt::xarray<double>;
  using Mat = xt::xarray<double>;
  using Arr = xt::xarray<double>;
  using shape_type = Vec::shape_type;

private:
  Arr& _F;

public:
  explicit lmi_oracle(Arr& F) : _F{F} {}

  auto chk_mtx(Mat &A, const Vec &x) const {
    // auto sub_mat = [](const Mat &M, size_t p) {
    //   return bnu::project(M, bnu::range(0, p), bnu::range(0, p));
    // };
    using xt::placeholders::_;
    using xt::linalg::dot;

    auto n = x.size();
    auto g = Vec(n);
    auto fj = -1.0;
    for (auto i = 0u; i < n; ++i) {
      auto Fi = xt::view(_F, xt::all(), xt::all(), i);
      A += Fi * x[i];
    }
    chol_ext Q(A);
    if (Q.is_sd()) {
      return std::make_tuple(g, fj);
    }
    auto v = Q.witness();
    auto p = v.size();
    auto sub = xt::range(_, p);
    fj = -dot(v, dot(xt::view(A, sub, sub), v))();
    for (auto i = 0u; i < n; ++i) {
      auto Fi = xt::view(_F, xt::all(), xt::all(), i);
      g[i] = -dot(v, dot(xt::view(Fi, sub, sub), v))();
    }
    return std::make_tuple(g, fj);
  }

  auto chk_spd_t(const Vec &x, double t) const {
    Mat A = _F[x.size()];
    auto m = A.size();
    for (auto i = 0u; i < m; ++i) {
      A(i, i) += t;
    }
    return this->chk_mtx(A, x);
  }

  auto chk_spd(const Vec &x) const {
    Mat A = _F[x.size()];
    return this->chk_mtx(A, x);
  }

  auto operator()(const Vec &x, double t) const {
    auto [g, fj] = this->chk_spd_t(x, t);
    if (fj < 0.0) {
      t -= 1.0;
    }
    return std::make_tuple(g, fj, t);
  }
};

#endif