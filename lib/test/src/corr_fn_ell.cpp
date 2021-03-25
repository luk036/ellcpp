// -*- coding: utf-8 -*-
#include <doctest/doctest.h>
#include <xtensor/xarray.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

extern Arr create_2d_sites(size_t, size_t);
extern Arr create_2d_isotropic(const Arr&, size_t);
extern std::tuple<Arr, size_t, bool> lsq_corr_poly2(
    const Arr&, const Arr&, size_t);
extern std::tuple<Arr, size_t, bool> mle_corr_poly(
    const Arr&, const Arr&, size_t);

TEST_CASE("check create_2d_isotropic")
{
    const auto s = create_2d_sites(5, 4);
    const auto Y = create_2d_isotropic(s, 3000);
    CHECK(s(6, 0) == doctest::Approx(2.5));
}

TEST_CASE("lsq_corr_fn")
{
    const auto s = create_2d_sites(10, 8);
    const auto Y = create_2d_isotropic(s, 3000);
    const auto [a, num_iters, feasible] = lsq_corr_poly2(Y, s, 4);
    CHECK(a[0] >= 0.);
    CHECK(feasible);
    CHECK(num_iters >= 679);
    CHECK(num_iters <= 712);
}

TEST_CASE("mle_corr_fn")
{
    const auto s = create_2d_sites(10, 8);
    const auto Y = create_2d_isotropic(s, 3000);
    const auto [a, num_iters, feasible] = mle_corr_poly(Y, s, 4);
    CHECK(a[0] >= 0.);
    CHECK(feasible);
    CHECK(num_iters >= 149);
    CHECK(num_iters <= 248);
}
