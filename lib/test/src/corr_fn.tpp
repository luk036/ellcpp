// -*- coding: utf-8 -*-
#include <catch.hpp>
#include <iostream>
#include <tuple>
#include <cutting_plane.hpp>
#include <ell.hpp>
#include <qmi_oracle.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xmath.hpp>
#include <vector>


// a fake dataset to make the bumps with
static const auto nx = 10;   // number of points
static const auto ny = 8;
static const auto n = nx*ny;
static const double s_end[] = {10., 8.};
static const auto sdkern = 0.1;  // width of kernel
static const auto var = 2.;     // standard derivation
static const auto tau = 0.00001;    // standard derivation of white noise
static const auto N = 200;  // number of samples


TEST_CASE("Corr_fn", "[corr_fn]") {
    using xt::linalg::dot;
    using Arr = xt::xarray<double>;

    xt::random::seed(5);
    // create sites s
    auto sx = xt::linspace(0, s_end[0], nx);
    auto sy = xt::linspace(0, s_end[1], ny);
    auto [xx, yy] = xt::meshgrid(sx, sy);
    // s = zip(xx.flatten(), yy.flatten())
    auto s = xt::vstack([xx.flatten(), yy.flatten()]).T;

    auto Sig = xt::ones<double>({n, n});
    for (auto i=0u; i<n; ++i) {
        for (auto j=i; j<n; ++j) {
            auto d = xt::array(s[j]) - xt::array(s[i]);
            Sig[i, j] = xt::exp(-sdkern * dot(d, d));
            Sig[j, i] = Sig[i, j];
        }
    }

    Arr A = xt::linalg::cholesky(Sig);
    Arr Ys = xt::zeros<double>({n, N});


    // ym = xt::random.randn(n)
    for (auto k=0u; k<N; ++k) {
        auto x = var * xt::random.randn(n);
        auto y = dot(A, x) + tau*xt::random.randn(n);
        Ys[:, k] = y;
    }

    Y = xt::cov(Ys, bias=True);

    // lsq_corr_bspline(Y, s, 4);
    lsq_corr_poly(Y, s, 4);
    
}
