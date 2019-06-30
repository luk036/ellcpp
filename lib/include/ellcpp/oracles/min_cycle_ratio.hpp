// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_MIN_CYCLE_RATIO_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_MIN_CYCLE_RATIO_HPP 1

#include <algorithm>
#include <py2cpp/py2cpp.hpp>
#include <tuple>
#include "parametric.hpp" // import max_parametric

// #include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

// template <typename Graph, typename Dict, typename T>
// void set_default(const Graph &G, Dict &map, T value) {
//     for (const auto& e : G.edges()) {
//         if (!map.contains(&e)) {
//             map[&e] = value;
//         }
//     }
// }

/*!
 * @brief
 *
 * @tparam Graph
 * @tparam Fn1
 * @tparam Fn2
 * @tparam T
 * @param G
 * @param get_cost
 * @param get_time
 * @return auto
 */
template<typename Graph, typename Fn1, typename Fn2, typename T>
auto min_cycle_ratio(Graph& G, Fn1 get_cost, Fn2 get_time, T&& /*! dummy */)
{
    using edge_t = decltype(*std::begin(G.edges()));
    edge_t e0    = *std::begin(G.edges());

    using cost_t = decltype(get_cost(G, e0));
    using time_t = decltype(get_time(G, e0));

    cost_t c0 = get_cost(G, e0);
    time_t t0 = get_time(G, e0);

    auto calc_weight = [&](const Graph&, T r, const auto& e) {
        return get_cost(G, e) - r * get_time(G, e);
    };

    auto calc_ratio = [&](const Graph& G, auto& C) {
        cost_t total_cost = cost_t(0);
        time_t total_time = time_t(0);
        for (const auto& e : C)
            total_cost += get_cost(G, e);
        for (const auto& e : C)
            total_time += get_time(G, e);
        return T(total_cost) / total_time;
    };

    // auto max_cost = *std::max_element(cost.begin(), cost.end());
    // auto min_time = *std::min_element(time.begin(), time.end());
    cost_t max_cost = c0;
    time_t min_time = t0;

    for (const auto& e : G.edges())
    {
        cost_t c = get_cost(G, e);
        time_t t = get_time(G, e);
        // std::cout << "mincost: c = " << c << '\n';
        if (max_cost < c) max_cost = c;
        if (min_time > t) min_time = t;
    }
    const auto r0 = T(max_cost * G.number_of_edges()) / min_time;
    return max_parametric(G, r0, calc_weight, calc_ratio);
}

#endif
