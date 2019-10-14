// -*- coding: utf-8 -*-
#pragma once

#include "parametric.hpp" // import max_parametric
#include <algorithm>
#include <numeric>
#include <py2cpp/py2cpp.hpp>

/*!
 * @brief
 *
 * @tparam Graph directed graph
 * @tparam Fn1
 * @tparam Fn2
 * @tparam T
 * @param G
 * @param get_cost
 * @param get_time
 * @return auto
 */
template <typename Graph, typename Fn1, typename Fn2, typename Container>
auto min_cycle_ratio(Graph& G, Fn1 get_cost, Fn2 get_time, Container& dist)
{
    using edge_t = typename Graph::edge_t;
    using T = typename Container::value_type;

    edge_t e0;
    for (auto e : G.edges())
    {
        e0 = e; // get the first edge (better idea???)
        break;
    }

    // auto max_cost = *std::max_element(cost.begin(), cost.end());
    // auto min_time = *std::min_element(time.begin(), time.end());
    auto max_cost = get_cost(G, e0);
    auto min_time = get_time(G, e0);
    for (auto e : G.edges())
    {
        auto c = get_cost(G, e);
        auto t = get_time(G, e);
        if (max_cost < c)
            max_cost = c;
        if (min_time > t)
            min_time = t;
    }
    auto r0 = T(max_cost * G.number_of_edges()) / min_time;

    using cost_T = decltype(get_cost(G, e0));
    using time_T = decltype(get_time(G, e0));

    auto calc_ratio = [&](const Graph& G, auto& C) {
        auto total_cost = cost_T(0);
        auto total_time = time_T(0);
        for (auto e : C)
        {
            total_cost += get_cost(G, e);
            total_time += get_time(G, e);
        }
        return T(total_cost) / total_time;
    };

    auto calc_weight = [&](const Graph&, T r, const auto& e) {
        return get_cost(G, e) - r * get_time(G, e);
    };
    return max_parametric(G, r0, calc_weight, calc_ratio, dist);
}
