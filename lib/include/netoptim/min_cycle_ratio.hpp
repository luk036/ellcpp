// -*- coding: utf-8 -*-
#pragma once

#include "parametric.hpp" // import max_parametric
#include <algorithm>
#include <numeric>
#include <py2cpp/py2cpp.hpp>

/*!
 * @brief minimum cost-to-time cycle ratio problem
 *
 *    This function solves the following network parametric problem:
 *
 *        max  r
 *        s.t. dist[v] − dist[u] ≤ cost(u, v) − r * time(u, v)
 *             ∀ e(u, v) ∈ G(V, E)
 *
 * @tparam Graph
 * @tparam Fn1
 * @tparam Fn2
 * @tparam Container
 * @param G
 * @param get_cost
 * @param get_time
 * @param dist
 * @return auto
 */
template <typename Graph, typename T, typename Fn1, typename Fn2,
    typename Container>
auto min_cycle_ratio(const Graph& G, T& r0, Fn1&& get_cost, Fn2&& get_time,
    Container&& dist, size_t max_iter = 1000)
{
    using edge_t = typename Graph::edge_t;
    using cost_T = decltype(get_cost(std::declval<edge_t>()));
    using time_T = decltype(get_time(std::declval<edge_t>()));

    auto calc_ratio = [&](const auto& C) -> T {
        auto total_cost = cost_T(0);
        auto total_time = time_T(0);
        for (auto&& e : C)
        {
            total_cost += get_cost(e);
            total_time += get_time(e);
        }
        return T(total_cost) / total_time;
    };

    auto calc_weight = [&](const T& r, const auto& e) {
        return get_cost(e) - r * get_time(e);
    };

    return max_parametric(G, r0, std::move(calc_weight), std::move(calc_ratio),
        std::forward<Container>(dist), max_iter);
}
