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

// extern std::size_t lsq_corr_poly(const Arr &, const Arr &, std::size_t);
extern std::tuple<size_t, bool> lsq_corr_poly2(const Arr &, const Arr &, std::size_t);
extern std::tuple<Arr, Arr> create_2d_isotropic(size_t, size_t, size_t);
extern std::tuple<size_t, bool> mle_corr_poly(const Arr &, const Arr &, std::size_t);
extern Arr construct_distance_matrix(const Arr &);

TEST_CASE("check create_2d_isotropic", "[corr_fn]") {
    auto [Y, s] = create_2d_isotropic(5, 4, 3000);
    // CHECK(Y(2,3) == Approx(1.9365965488224368));
    CHECK(s(6, 0) == Approx(2.5));
    auto D1 = construct_distance_matrix(s);
    CHECK(D1(2, 4) == Approx(5.0));    
}

TEST_CASE("lsq_corr_fn", "[corr_fn]") {
    auto [Y, s] = create_2d_isotropic(5, 4, 3000);
    // lsq_corr_bspline(Y, s, 4);
    auto [num_iters, feasible] = lsq_corr_poly2(Y, s, 4);
    CHECK(feasible);
    CHECK(num_iters >= 8);
    CHECK(num_iters <= 657);
}

TEST_CASE("mle_corr_fn", "[corr_fn]") {
    auto [Y, s] = create_2d_isotropic(5, 4, 3000);
    // lsq_corr_bspline(Y, s, 4);
    auto [num_iters, feasible] = mle_corr_poly(Y, s, 4);
    CHECK(feasible);
    CHECK(num_iters >= 50);
    CHECK(num_iters <= 657);
}
