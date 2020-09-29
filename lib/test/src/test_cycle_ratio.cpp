// -*- coding: utf-8 -*-
#include <array>
#include <doctest.h>
#include <netoptim/min_cycle_ratio.hpp>
#include <py2cpp/fractions.hpp> // import Fraction
#include <xnetwork/classes/digraphs.hpp>


TEST_CASE("Test Cycle Ratio")
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

    auto dist = std::vector<fun::Fraction<int>>(
        G.number_of_nodes(), fun::Fraction<int>(0));
    auto r = fun::Fraction<int>(5);
    const auto c = min_cycle_ratio(G, r, get_cost, get_time, dist);
    CHECK(!c.empty());
    CHECK(c.size() == 5);
    CHECK(r == fun::Fraction<int>(9, 5));
}

TEST_CASE("Test Cycle Ratio of Timing Graph")
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

    auto dist = std::vector<fun::Fraction<int>>(
        G.number_of_nodes(), fun::Fraction<int>(0));
    auto r = fun::Fraction<int>(7);
    const auto c = min_cycle_ratio(G, r, get_cost, get_time, dist);
    CHECK(!c.empty());
    CHECK(r == fun::Fraction<int>(1, 1));
    CHECK(c.size() == 3);
    // print(r);
    // print(c);
    // print(dist.items());
}
