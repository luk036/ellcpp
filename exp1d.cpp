#include "cutting_plane.hpp"
#include "ell.hpp"
#include <xtensor/xarray.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <iostream>
#include <cmath>

int main() {
  using xt::linalg::dot;
  using Arr = xt::xarray<double, xt::layout_type::row_major>;

  const auto n = 4u;   // number of points
  const auto s_begin = 1U;
  const auto s_end = 10u;
  const auto sdkern = 0.5;  // width of kernel
  const auto var = 2.;     // standard derivation
  const auto N = 500u;      // number of samples

  // auto dist = s_end - s_begin;
  Arr s = xt::linspace(s_begin, s_end, n);

  //shape_type shape = {n, n};
  Arr Sig = xt::ones<double>({n, n});

  for (auto i = 0U; i < n; ++i) {
    for (auto j = i+1; j < n; ++j) {
      const auto d = s[j] - s[i];
      Sig(i, j) = std::exp(-0.5*(d*d)/(sdkern*sdkern)/2);
      Sig(j, i) = Sig(i,j);
    }
  }

  Arr A = xt::linalg::cholesky(Sig);
  Arr Ys = zeros({n, N});

  Arr ym = xt::random::randn<double>({n});
  for (auto k = 0U; k < N; ++k) {
    Arr x = var * xt::random::randn<double>({n});
    Arr y = dot(A,x) + ym + 0.5*xt::random::randn<double>({n});
    xt::view(Ys, xt::all(), k) = y;
  }
  
  Arr Ys_T = xt::transpose(Ys);
  Arr Y = dot(Ys, Ys_T);
  Y *= 1./N;

  //auto Y = xt::cov(Ys, bias=True);
  // std::cout << Y << std::endl;
}
