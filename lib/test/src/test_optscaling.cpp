// -*- coding: utf-8 -*-
#include <algorithm>
#include <catch2/catch.hpp>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/ell1d.hpp>
#include <ellcpp/oracles/optscaling3_oracle.hpp> // import optscaling
#include <ellcpp/oracles/optscaling_oracle.hpp>  // import optscaling
#include <xnetwork/classes/digraphs.hpp>
#include <utility> // for std::pair
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
    const auto num_of_nodes = sizeof(elem) / sizeof(double);
    double cost[num_of_nodes];
    for (size_t i = 0U; i < num_of_nodes; ++i)
    {
        cost[i] = std::log(std::abs(elem[i]));
    }

    auto get_cost = [&](const auto& G, const auto& e) -> int {
        auto [u, v] = G.end_points(e);
        return cost[G[u][v]];
    };


    auto [cmin, cmax] = std::minmax_element(cost, cost + num_of_nodes);
    auto x0 = Arr {*cmax, *cmin};
    auto t = *cmax - *cmin;
    auto E = ell {1.5 * t, x0};
    auto dist = std::vector(G.number_of_nodes(), 0.);

    auto P = optscaling_oracle {G, dist, get_cost};
    auto [_, ell_info] =
        cutting_plane_dc(P, E, std::numeric_limits<double>::max());

    CHECK(ell_info.feasible);
    CHECK(ell_info.num_iters <= 27);
}

TEST_CASE("Test Optimal Scaling (binary search)", "[test_optscaling]")
{
    auto G = create_test_case1();

    const auto elem = std::array {1.2, 2.3, 3.4, -4.5, 5.6};
    const auto num_of_nodes = sizeof(elem) / sizeof(double);
    double cost[num_of_nodes];
    for (size_t i = 0U; i < num_of_nodes; ++i)
    {
        cost[i] = std::log(std::abs(elem[i]));
    }

    auto get_cost = [&](const auto& G, const auto& e) -> int {
        auto [u, v] = G.end_points(e);
        return cost[G[u][v]];
    };


    auto [cmin, cmax] = std::minmax_element(cost, cost + num_of_nodes);

    auto Iv = ell1d(*cmin, *cmax);
    auto dist = std::vector(G.number_of_nodes(), 0.);

    auto Q = optscaling3_oracle {G, dist, get_cost};
    auto P = bsearch_adaptor{Q, Iv};
    auto bs_info = bsearch(P, std::tuple {0., 1.001 * (*cmax - *cmin)});

    CHECK(bs_info.feasible);
    CHECK(bs_info.num_iters <= 27);
}
