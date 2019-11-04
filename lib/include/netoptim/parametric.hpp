// -*- coding: utf-8 -*-
#pragma once

#include "neg_cycle.hpp" // import negCycleFinder
#include <iostream>
#include <tuple>
#include <vector>

/*!
 * @brief maximum parametric problem
 *
 *    max  r
 *    s.t. dist[v] - dist[u] <= d(u, v, r)
 *         for all (u, v) : G
 *
 * @tparam Graph
 * @tparam T
 * @tparam Fn1
 * @tparam Fn2
 * @tparam Container
 * @param G directed graph
 * @param[inout] r_opt parameter to be maximized, initially a large number
 * @param d monotone decreasing function w.r.t. r
 * @param zero_cancel
 * @param[inout] dist
 * @return optimal r and the critical cycle
 */
template <typename Graph, typename T, typename Fn1, typename Fn2,
    typename Container>
auto max_parametric(
    Graph& G, T r_opt, Fn1& d, Fn2& zero_cancel, Container& dist)
{
    using edge_t = typename Graph::edge_t;

    auto get_weight = [&](const edge_t& e) -> T { // int???
        return d(G, r_opt, e);
    };

    auto S = negCycleFinder(G);
    auto C_opt = std::vector<edge_t> {}; // should initial outside

    while (true)
    {
        const auto& C_min = S.find_neg_cycle(dist, get_weight);
        if (C_min.empty())
        {
            break;
        }

        const auto& r_min = zero_cancel(G, C_min);
        if (r_min >= r_opt)
        {
            break;
        }

        C_opt = C_min;
        r_opt = r_min;

        // update ???
        for (const auto& e : C_opt)
        {
            auto [u, v] = G.end_points(e);
            dist[u] = dist[v] - get_weight(e);
        }
    }

    return std::tuple {r_opt, std::move(C_opt)};
}
