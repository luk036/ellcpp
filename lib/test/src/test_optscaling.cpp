// -*- coding: utf-8 -*-
#include <algorithm>
#include <array>
#include <doctest.h>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
// #include <ellcpp/ell1d.hpp>
#include <cstddef>
#include <ellcpp/oracles/optscaling_oracle.hpp> // import optscaling
#include <limits>
#include <utility> // for std::pair
#include <xnetwork/classes/digraphs.hpp>
#include <xtensor/xarray.hpp>

TEST_CASE("Test Optimal Scaling (two varaibles)")
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
        std::array <Edge, 5>{Edge(A, B), Edge(B, C), Edge(C, D), Edge(D, E), Edge(E, A)};
    const auto indices = std::array<int, 5> {0, 1, 2, 3, 4};
    auto G = xn::SimpleDiGraphS(py::range<int>(num_nodes));
    G.add_edges_from(edges, indices);

    const auto elem = std::array<double, 5> {1.2, 2.3, 3.4, -4.5, 5.6};
    constexpr auto num_of_nodes = sizeof(elem) / sizeof(double);

    // @todo: use std::array<> instead
    // double cost[num_of_nodes];
    auto cost = std::array<double, num_of_nodes> {};
    for (size_t i = 0U; i != num_of_nodes; ++i)
    {
        cost[i] = std::log(std::abs(elem[i]));
    }

    auto get_cost = [&](const auto& edge) -> double {
        const auto e = G.end_points(edge);
        return cost[G[e.first].at(e.second)];
    };


    const auto [cmin, cmax] =
        std::minmax_element(std::begin(cost), std::end(cost));
    auto x0 = Arr {*cmax, *cmin};
    auto t1 = *cmax - *cmin;
    auto El = ell {1.5 * t1, std::move(x0)};
    auto dist = std::vector<double>(G.number_of_nodes(), 0.);
    auto P = optscaling_oracle<xn::SimpleDiGraphS, std::vector<double>, 
    decltype(get_cost)> {G, dist, get_cost};
    auto t = 1.e100; // std::numeric_limits<double>::max()
    const auto [x, ell_info] = cutting_plane_dc(std::move(P), std::move(El), t);

    CHECK(x[1] >= x[0]);
    CHECK(ell_info.feasible);
    CHECK(ell_info.num_iters <= 27);
}
