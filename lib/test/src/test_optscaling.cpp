// -*- coding: utf-8 -*-
#include <algorithm>
#include <array>
#include <catch2/catch.hpp>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
// #include <ellcpp/ell1d.hpp>
#include <ellcpp/oracles/optscaling_oracle.hpp> // import optscaling
#include <utility>                              // for std::pair
#include <xnetwork/classes/digraphs.hpp>
#include <xtensor/xarray.hpp>

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
    const auto edges =
        std::array {Edge(A, B), Edge(B, C), Edge(C, D), Edge(D, E), Edge(E, A)};
    const auto indices = std::array {0, 1, 2, 3, 4};
    auto g = xn::DiGraphS(py::range<int>(num_nodes));
    g.add_edges_from(edges, indices);
    return g;
}

TEST_CASE("Test Optimal Scaling (two varaibles)", "[test_optscaling]")
{
    auto G = create_test_case1();

    const auto elem = std::array {1.2, 2.3, 3.4, -4.5, 5.6};
    constexpr auto num_of_nodes = sizeof(elem) / sizeof(double);

    // @todo: use std::array<> instead
    // double cost[num_of_nodes];
    auto cost = std::array<double, num_of_nodes> {};
    for (size_t i = 0U; i < num_of_nodes; ++i)
    {
        cost[i] = std::log(std::abs(elem[i]));
    }

    auto get_cost = [&](const auto& e) -> int {
        const auto [u, v] = G.end_points(e);
        return cost[G[u][v]];
    };


    const auto [cmin, cmax] =
        std::minmax_element(std::begin(cost), std::end(cost));
    const auto x0 = Arr {*cmax, *cmin};
    auto t1 = *cmax - *cmin;
    auto E = ell {1.5 * t1, x0};
    auto dist = std::vector(G.number_of_nodes(), 0.);

    auto P = optscaling_oracle {G, dist, get_cost};
    auto t = std::numeric_limits<double>::max();
    [[maybe_unused]] const auto [_, ell_info] = cutting_plane_dc(P, E, t);

    CHECK(ell_info.feasible);
    CHECK(ell_info.num_iters <= 27);
}
