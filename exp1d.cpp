#include "cutting_plane.hpp"
#include "ell.hpp"
#include "profit_oracle.hpp"
#include <xtensor/xarray.hpp>
#include <xtensor/xrandom.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <iostream>
#include <cmath>

int main() {
  using xt::linalg::dot;

  auto n = 4u;   // number of points
  auto s_begin = 1u;
  auto s_end = 10u;
  auto sdkern = 0.5;  // width of kernel
  auto var = 2.0;     // standard derivation
  auto N = 500u;      // number of samples

  auto dist = s_end - s_begin;
  auto s = xt::linspace(s_begin, s_end, n);

  //shape_type shape = {n, n};
  xt::xarray<double> Sig = xt::ones<double>({n, n});

  for (auto i = 0u; i < n; ++i) {
    for (auto j = i+1; j < n; ++j) {
      auto d = s[j] - s[i];
      Sig(i, j) = std::exp(-0.5*(d*d)/(sdkern*sdkern)/2.0);
      Sig(j, i) = Sig(i,j);
    }
  }

  auto A = xt::linalg::cholesky(Sig);
  xt::xarray<double> Ys = xt::zeros<double>({n, N});

  auto ym = xt::random::randn<double>({n});
  for (auto k = 0u; k < N; ++k) {
    auto x = var * xt::random::randn<double>({n});
    auto y = dot(A,x) + ym + 0.5*xt::random::randn<double>({n});
    xt::view(Ys, xt::all(), k) = y;
  }
  
  auto Ys_T = xt::transpose(Ys);
  auto Y = dot(Ys, Ys_T);
  Y *= 1./N;

  //auto Y = xt::cov(Ys, bias=True);
  std::cout << Y << std::endl;

}