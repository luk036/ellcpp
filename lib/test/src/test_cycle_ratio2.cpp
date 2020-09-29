// -*- coding: utf-8 -*-
#include <array>
#include <cstddef>
#include <doctest.h>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell1d.hpp>
#include <ellcpp/oracles/cycle_ratio_oracle.hpp> // import cycle_ratio
#include <limits>
#include <utility> // for std::pair
#include <xnetwork/classes/digraphs.hpp>

TEST_CASE("Test Cycle Ratio 2")
{
    using Edge = std::pair<int, int>;
    constexpr auto num_nodes = 5;
    enum nodes
    {
        A,
        B,
        C,
        D,
        E
    };
    const auto edges = std::array<Edge, 5> {
        Edge {A, B}, Edge {B, C}, Edge {C, D}, Edge {D, E}, Edge {E, A}};
    const auto indices = std::array<int, 5> {0, 1, 2, 3, 4};
    auto G = xn::SimpleDiGraphS {py::range<int>(num_nodes)};
    G.add_edges_from(edges, indices);

    const auto cost = std::array<int, 5> {5, 1, 1, 1, 1};

    const auto get_cost = [&](const auto& edge) -> int {
        const auto e = G.end_points(edge);
        return cost[G[e.first].at(e.second)];
    };
    const auto get_time = [&](const auto & /*e*/) -> int { return 1; };

    auto dist = std::vector<double>(G.number_of_nodes(), 0.);
    auto Ell = ell1d {-100., 100.};
    auto P = cycle_ratio_oracle<xn::SimpleDiGraphS, std::vector<double>,
        decltype(get_cost), decltype(get_time)> {G, dist, get_cost, get_time};
    auto r = -1.e100; // std::numeric_limits<double>::min()
    const auto opts = Options {2000, 1e-12};
    const auto [x, ell_info] = cutting_plane_dc(P, Ell, r, opts);
    CHECK(ell_info.feasible);
    CHECK(r == doctest::Approx(9. / 5.));
}

TEST_CASE("Test Cycle Ratio of Timing Graph 2")
{
    using Edge = std::pair<int, int>;
    constexpr auto num_nodes = 3;
    enum nodes
    {
        A,
        B,
        C
    };
    const auto edges = std::array<Edge, 6> {Edge {A, B}, Edge {B, A},
        Edge {B, C}, Edge {C, B}, Edge {C, A}, Edge {A, C}};
    // make sure no parallel edges!!!

    const auto indices = std::array<int, 6> {0, 1, 2, 3, 4, 5};
    auto G = xn::SimpleDiGraphS {py::range<int>(num_nodes)};
    G.add_edges_from(edges, indices);

    const auto cost = std::array<int, 6> {7, -1, 3, 0, 2, 4};

    const auto get_cost = [&](const auto& edge) -> int {
        const auto e = G.end_points(edge);
        return cost[G[e.first].at(e.second)];
    };
    const auto get_time = [&](const auto & /*e*/) -> int { return 1; };

    auto dist = std::vector<double>(G.number_of_nodes(), 0.);
    auto E = ell1d {-100., 100.};
    auto P = cycle_ratio_oracle<xn::SimpleDiGraphS, std::vector<double>,
        decltype(get_cost), decltype(get_time)> {G, dist, get_cost, get_time};
    auto r = -1.e100; // std::numeric_limits<double>::min()
    const auto opts = Options {2000, 1e-12};
    const auto [x, ell_info] = cutting_plane_dc(P, E, r, opts);
    CHECK(ell_info.feasible);
    CHECK(r == doctest::Approx(1.));
    // print(r);
    // print(c);
    // print(dist.items());
}
