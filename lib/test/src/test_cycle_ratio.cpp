// -*- coding: utf-8 -*-
#include <array>
#include <catch2/catch.hpp>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell1d.hpp>
#include <ellcpp/oracles/cycle_ratio_oracle.hpp> // import cycle_ratio
#include <utility>                              // for std::pair
#include <xnetwork/classes/digraphs.hpp>

/*!
 * @brief Create a test case1 object
 *
 * @return auto
 */
static auto create_test_case1()
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
    const auto edges = std::array {
        Edge {A, B}, Edge {B, C}, Edge {C, D}, Edge {D, E}, Edge {E, A}};
    const auto indices = std::array {0, 1, 2, 3, 4};
    auto g = xn::DiGraphS {py::range<int>(num_nodes)};
    g.add_edges_from(edges, indices);
    return g;
}

/*!
 * @brief Create a test case timing object
 *
 * @return auto
 */
static auto create_test_case_timing()
{
    using Edge = std::pair<int, int>;
    constexpr auto num_nodes = 3;
    enum nodes
    {
        A,
        B,
        C
    };
    const auto edges = std::array {Edge {A, B}, Edge {B, A}, Edge {B, C},
        Edge {C, B}, Edge {C, A}, Edge {A, C}};
    // make sure no parallel edges!!!

    const auto indices = std::array {0, 1, 2, 3, 4, 5};
    auto g = xn::DiGraphS {py::range<int>(num_nodes)};
    g.add_edges_from(edges, indices);
    return g;
}

TEST_CASE("Test Cycle Ratio", "[test_cycle_ratio]")
{
    const auto G = create_test_case1();
    const auto cost = std::array {5, 1, 1, 1, 1};

    const auto get_cost = [&](const auto& e) -> int {
        auto [u, v] = G.end_points(e);
        return cost[G[u][v]];
    };
    const auto get_time = [&](const auto & /*e*/) -> int { return 1; };

    auto dist = std::vector(G.number_of_nodes(), 0.);
    auto E = ell1d {-100., 100.};
    auto P = cycle_ratio_oracle {G, dist, get_cost, get_time};
    auto r = std::numeric_limits<double>::min();
    const auto opts = Options{2000, 1e-12};
    const auto [x, ell_info] = cutting_plane_dc(P, E, r, opts);
    CHECK(ell_info.feasible);
    CHECK(r == Approx(9. / 5.));
}

TEST_CASE("Test Cycle Ratio of Timing Graph", "[test_cycle_ratio]")
{
    const auto G = create_test_case_timing();
    const auto cost = std::array {7, -1, 3, 0, 2, 4};

    const auto get_cost = [&](const auto& e) -> int {
        auto [u, v] = G.end_points(e);
        return cost[G[u][v]];
    };
    const auto get_time = [&](const auto & /*e*/) -> int { return 1; };

    auto dist = std::vector(G.number_of_nodes(), 0.);
    auto E = ell1d {-100., 100.};
    auto P = cycle_ratio_oracle {G, dist, get_cost, get_time};
    auto r = std::numeric_limits<double>::min();
    const auto opts = Options{2000, 1e-12};
    const auto [x, ell_info] = cutting_plane_dc(P, E, r, opts);
    CHECK(ell_info.feasible);
    CHECK(r == Approx(1.));
    // print(r);
    // print(c);
    // print(dist.items());
}
