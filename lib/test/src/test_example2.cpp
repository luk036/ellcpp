/* -*- coding: utf-8 -*- */
#include <catch.hpp>
// #include <iostream>
// #include <tuple>

#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

auto my_oracle2(const Arr& z)
{
    auto x = z[0];
    auto y = z[1];

    // constraint 1: x + y <= 3
    auto fj = x + y - 3.;
    if (fj > 0) { return std::tuple{Arr{1., 1.}, fj, false}; }

    // constraint 2: x - y >= 1
    fj = -x + y + 1.;
    if (fj > 0) { return std::tuple{Arr{-1., 1.}, fj, false}; }

    return std::tuple{Arr{0., 0.}, 0., true};
}

TEST_CASE("Example 2", "[example2]")
{
    auto x0 = Arr{0., 0.}; // initial x0
    auto E  = ell{10., x0};
    auto P  = my_oracle2;

    auto ell_info = cutting_plane_feas(P, E);
    CHECK(ell_info.feasible);
    // std::cout << "Example 2 result: " << niter << "," << feasible << "," <<
    // status << "\n"; std::cout << "Example 2 xbest: " << xb << "\n";
}