// -*- coding: utf-8 -*-
#include <catch.hpp>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <iostream>
// #include <ellcpp/oracles/qmi_oracle.hpp>
#include <tuple>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xrandom.hpp>

using Arr = xt::xarray<double>;

extern std::size_t lsq_corr_poly(const Arr &, const Arr &, std::size_t);
extern std::size_t lsq_corr_poly2(const Arr &, const Arr &, std::size_t);

// a fake dataset to make the bumps with
static const auto nx = 11; // number of points
static const auto ny = 9;
static const auto n = nx * ny;
static const double s_end[] = {10., 8.};
static const auto sdkern = 0.1;  // width of kernel
static const auto var = 2.;      // standard derivation
static const auto tau = 0.00001; // standard derivation of white noise
static const auto N = 200;       // number of samples

TEST_CASE("Corr_fn", "[corr_fn]") {
    using xt::linalg::dot;

    xt::random::seed(5);
    // create sites s
    Arr sx = xt::linspace<double>(0., 10., nx);
    Arr sy = xt::linspace<double>(0., 8., ny);
    auto [xx, yy] = xt::meshgrid(sx, sy);
    // s = zip(xx.flatten(), yy.flatten())
    Arr s = xt::stack(xt::xtuple(xt::flatten(xx), xt::flatten(yy)), 0);
    s = xt::transpose(s);

    Arr Sig = xt::zeros<double>({n, n});
    for (auto i = 0u; i < n; ++i) {
        for (auto j = i; j < n; ++j) {
            Arr d = xt::view(s, j, xt::all()) - xt::view(s, i, xt::all());
            double g = -sdkern * dot(d, d)();
            Sig(i, j) = std::exp(g);
            Sig(j, i) = Sig(i, j);
        }
    }

    Arr A = xt::linalg::cholesky(Sig);
    // Arr Ys = xt::zeros<double>({n, N});

    // ym = xt::random.randn(n)
    // Y = xt::cov(Ys, bias=True);
    Arr Y = xt::zeros<double>({n, n});
    for (auto k = 0u; k < N; ++k) {
        Arr x = var * xt::random::randn<double>({n});
        // auto y = dot(A, x)() + tau*xt::random::randn<double>({n});
        Arr y = dot(A, x);
        // Ys[:, k] = y;
        Y += xt::linalg::outer(y, y);
    }
    Y /= N;

    // lsq_corr_bspline(Y, s, 4);
    auto num_iters = lsq_corr_poly2(Y, s, 4);
    CHECK(num_iters >= 88);
    CHECK(num_iters <= 657);
}
