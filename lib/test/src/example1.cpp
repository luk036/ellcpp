// -*- coding: utf-8 -*-
#include <catch.hpp>
#include <iostream>
#include <tuple>

#include <cutting_plane.hpp>
#include <ell.hpp>

using Arr = xt::xarray<double>;

auto my_oracle(const Arr &z, double t) {
    auto x = z[0], y = z[1];

    // constraint 1: x + y <= 3
    auto fj = x + y - 3;
    if (fj > 0.) {
        return std::tuple{Arr{1., 1.}, fj, t};
    }
    // constraint 2: x - y >= 1
    fj = -x + y + 1;
    if (fj > 0.) {
        return std::tuple{Arr{-1., 1.}, fj, t};
    }
    // objective: maximize x + y
    auto f0 = x + y;
    fj = t - f0;
    if (fj < 0.) {
        fj = 0.;
        t = f0;
    }
    return std::tuple{Arr{-1., -1.}, fj, t};
}

TEST_CASE("Example 1", "[example1]") {
    auto x0 = Arr{0., 0.}; // initial x0
    auto E = ell(10., x0);
    auto P = my_oracle;
    auto [xb, fb, niter, feasible, status] = cutting_plane_dc(P, E, -100.);
    CHECK(feasible);
    std::cout << niter << "," << feasible << "," << status << "\n";
    std::cout << xb << "\n";
}