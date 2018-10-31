// -*- coding: utf-8 -*-
#include <catch.hpp>
#include <iostream>
#include <tuple>

#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>

using Arr = xt::xarray<double>;

class nocopymove {
  private:
    int i = 3;

  public:
    nocopymove() = default;
    nocopymove(const nocopymove &) = delete;
    nocopymove(nocopymove &&) = default;
};

auto func() {
    nocopymove a{};
    return std::tuple{std::move(a), 1.};
}

auto my_oracle(const Arr &z, double t) {
    double x = z[0], y = z[1];

    // constraint 1: x + y <= 3
    double fj = x + y - 3;
    if (fj > 0) {
        return std::tuple{Arr{1., 1.}, fj, t};
    }
    // constraint 2: x - y >= 1
    fj = -x + y + 1;
    if (fj > 0) {
        return std::tuple{Arr{-1., 1.}, fj, t};
    }
    // objective: maximize x + y
    double f0 = x + y;
    fj = t - f0;
    if (fj < 0) {
        fj = 0.;
        t = f0;
    }
    return std::tuple{Arr{-1., -1.}, fj, t};
}

TEST_CASE("Example 1", "[example1]") {
    Arr x0{0., 0.}; // initial x0
    ell E(10., x0);
    auto P = my_oracle;
    auto [xb, fb, niter, feasible, status] = cutting_plane_dc(P, E, -100.);
    CHECK(feasible);
    std::cout << "Example 1 result: " << fb << ", " << niter << "," << feasible << "," << status << "\n";
    std::cout << "Example 1 xbest: " << xb << "\n";

    auto [c, d] = func();
}
