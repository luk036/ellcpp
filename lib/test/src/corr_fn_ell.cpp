// -*- coding: utf-8 -*-
#include <catch.hpp>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <iostream>
// #include <ellcpp/oracles/qmi_oracle.hpp>
#include <tuple>
#include <vector>
#include <xtensor/xarray.hpp>

using Arr = xt::xarray<double>;

extern std::size_t lsq_corr_poly(const Arr &, const Arr &, std::size_t);
extern std::size_t lsq_corr_poly2(const Arr &, const Arr &, std::size_t);
extern std::tuple<Arr, Arr> create_2d_isotropic(size_t, size_t, size_t);
extern std::size_t mle_corr_poly(const Arr &, const Arr &, std::size_t);

TEST_CASE("lsq_corr_fn", "[corr_fn]") {
    auto [Y, s] = create_2d_isotropic(11, 9, 200);
    // lsq_corr_bspline(Y, s, 4);
    auto num_iters = lsq_corr_poly2(Y, s, 4);
    CHECK(num_iters >= 88);
    CHECK(num_iters <= 657);
}

TEST_CASE("mle_corr_fn", "[corr_fn]") {
    auto [Y, s] = create_2d_isotropic(5, 4, 3000);
    // lsq_corr_bspline(Y, s, 4);
    auto num_iters = lsq_corr_poly2(Y, s, 4);
    CHECK(num_iters >= 88);
    CHECK(num_iters <= 657);
}
