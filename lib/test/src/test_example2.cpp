/* -*- coding: utf-8 -*- */
#include <catch.hpp>
#include <iostream>
#include <tuple>

#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>

using Arr = xt::xarray<double>;

auto my_oracle2(const Arr &z) {
    double x = z[0], y = z[1];

    // constraint 1: x + y <= 3
    double fj = x + y - 3;
    if (fj > 0) {
        return std::tuple{Arr{1., 1.}, fj, false};
    }
    // constraint 2: x - y >= 1
    fj = -x + y + 1;
    if (fj > 0) {
        return std::tuple{Arr{-1., 1.}, fj, false};
    }
    return std::tuple{Arr{0., 0.}, 0., true};
}

TEST_CASE("Example 2", "[example2]") {
    Arr x0{0., 0.}; // initial x0
    ell E(10., x0);
    auto P = my_oracle2;
    auto [xb, niter, feasible, status] = cutting_plane_feas(P, E);
    CHECK(feasible);
    std::cout << "Example 2 result: " << niter << "," << feasible << "," << status << "\n";
    std::cout << "Example 2 xbest: " << xb << "\n";
}