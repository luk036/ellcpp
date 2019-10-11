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
 *    s.t. dist[v] - dist[v] <= d(u, v, r);
 *         for all (u, v) : G
 *
 * @tparam Graph
 * @tparam T
 * @tparam Fn1
 * @tparam Fn2
 * @param G directed graph
 * @param r parameter to be maximized, initially a large number
 * @param d monotone decreasing function w.r.t. r
 * @param zero_cancel
 * @return:
 *       r_opt -- optimal value
 *       C_opt -- Most critial cycle
 *       dist -- optimal sol"n
 */
template <typename Graph, typename T, typename Fn1, typename Fn2>
auto max_parametric(Graph& G, T r, Fn1& d, Fn2& zero_cancel)
{
    using edge_t = typename Graph::edge_t;

    auto get_weight = [d, r](const Graph& G,
                          const edge_t& e) -> T { // int???
        return d(G, r, e);
    };

    auto S = negCycleFinder(G, get_weight);
    auto C_opt = std::vector<edge_t> {};
    auto r_opt = r;

    while (true)
    {
        const auto& C_min = S.neg_cycle_relax();
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
            S._dist[u] = S._dist[v] - get_weight(G, e);
        }
    }

    return std::tuple {r_opt, std::move(C_opt)};
}
