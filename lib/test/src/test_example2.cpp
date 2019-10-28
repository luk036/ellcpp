/* -*- coding: utf-8 -*- */
#include <catch2/catch.hpp>
// #include <iostream>
// #include <tuple>

#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <optional>

using Arr = xt::xarray<double, xt::layout_type::row_major>;
using Cut = std::tuple<Arr, double>;

auto my_oracle2(const Arr& z) -> std::optional<Cut>
{
    auto x = z[0];
    auto y = z[1];

    // constraint 1: x + y <= 3
    auto fj = x + y - 3.;
    if (fj > 0.)
    {
        return {{Arr {1., 1.}, fj}};
    }

    // constraint 2: x - y >= 1
    fj = -x + y + 1.;
    if (fj > 0.)
    {
        return {{Arr {-1., 1.}, fj}};
    }

    return {};
}

TEST_CASE("Example 2", "[example2]")
{
    auto E = ell {10., Arr {0., 0.}};
    auto P = my_oracle2;

    auto ell_info = cutting_plane_feas(P, E);
    CHECK(ell_info.feasible);
    // std::cout << "Example 2 result: " << niter << "," << feasible << "," <<
    // status << "\n"; std::cout << "Example 2 xbest: " << xb << "\n";
}