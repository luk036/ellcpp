// -*- coding: utf-8 -*-
#pragma once

#include "neg_cycle.hpp" // import negCycleFinder
#include <iostream>
#include <tuple>
#include <vector>

/*!
 * @brief maximum parametric problem
 *
 *    This function solves the following network parametric problem:
 *
 *        max  r
 *        s.t. dist[v] − dist[u] ≤ d(u, v, r)
 *             ∀ e(u, v) ∈ G(V, E)
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
auto max_parametric(const Graph& G, T& r_opt, Fn1&& d, Fn2&& zero_cancel,
    Container&& dist, size_t max_iter = 1000)
{
    using edge_t = typename Graph::edge_t;

    auto get_weight = [&](const edge_t& e) -> T { // int???
        return d(r_opt, e);
    };

    auto S = negCycleFinder(G);
    auto C_opt = std::vector<edge_t> {}; // should initial outside

    auto niter = 0U;
    for (; niter < max_iter; ++niter)
    {
        const auto& C_min = S.find_neg_cycle(std::forward<Container>(dist),
                                             std::move(get_weight));
        if (C_min.empty())
        {
            break;
        }

        const auto& r_min = zero_cancel(C_min);
        if (r_min >= r_opt)
        {
            break;
        }

        C_opt = C_min;
        r_opt = r_min;

        // update ???
        for (auto&& e : C_opt)
        {
            const auto [u, v] = G.end_points(e);
            dist[u] = dist[v] - get_weight(e);
        }
    }

    return C_opt;
}
