#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI_ORACLE_HPP 1

//#include "mat.hpp"
#include "chol_ext.hpp"
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

inline static auto quad(const xt::xarray<double> &A,
                        const xt::xarray<double> &v, size_t p) {
  double res = 0.0;
  for (auto i = 0u; i < p; ++i) {
    for (auto j = 0u; j < p; ++j) {
      res += A(i, j) * v(i) * v(j);
    }
  }
  return res;
}

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
  Arr &_F;
  Arr &_F0;

public:
  explicit lmi_oracle(Arr &F, Arr &B) : _F{F}, _F0{B} {}

  auto chk_mtx(Arr A, const Arr &x) const {
    using xt::linalg::dot;
    using xt::placeholders::_;

    auto n = x.size();
    Arr g = xt::zeros<double>({n});
    auto fj = -1.0;
    for (auto i = 0u; i < n; ++i) {
      auto Fi = xt::view(_F, i, xt::all(), xt::all());
      // Arr Fi = _F(i);
      A -= Fi * x(i);
    }
    chol_ext Q(A);
    if (Q.is_sd()) {
      return std::tuple{g, fj};
    }
    Arr v = Q.witness();
    double p = v.size();
    fj = 1.0;
    for (auto i = 0u; i < n; ++i) {
      auto Fipp = xt::view(_F, i, xt::range(_, p), xt::range(_, p));
      g(i) = quad(Fipp, v, p);
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

  auto chk_spd(const Arr &x) const { return this->chk_mtx(_F0, x); }

  auto operator()(const Arr &x, double t) const {
    auto [g, fj] = this->chk_spd_t(x, t);
    if (fj < 0.0) {
      t -= 1.0;
    }
    return std::make_tuple(g, fj, t);
  }
};

#endif