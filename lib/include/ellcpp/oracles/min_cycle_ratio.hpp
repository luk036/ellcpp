// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_MIN_CYCLE_RATIO_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_MIN_CYCLE_RATIO_HPP 1

#include "parametric.hpp" // import max_parametric
#include <py2cpp/py2cpp.hpp>
#include <algorithm>

template <typename Graph, typename Dict, typename T>
void set_default(const Graph &G, Dict& map, T value) {
    for (auto e : G.edges()) {
        if (!map.contains(&e)) {
            map[&e] = value;
        }
    }
}


template <typename Graph, typename Dict>
auto min_cycle_ratio(const Graph &G, Dict& cost, Dict& time) {
    using edge_t = decltype(*(G.edges().begin()));

    auto calc_weight = [cost, time](const Graph &G, double r, edge_t& e) {
        return cost[&e] - r * time[&e];
    };

    auto calc_ratio = [cost, time](const Graph &G, auto& C) {
        auto total_cost = 0;
        auto total_time = 0;
        for (auto e : C) {
            total_cost += cost[e];
            total_time += time[e];
        }
        return total_cost / total_time;
    };

    set_default(G, cost, 1);
    set_default(G, time, 1);
    auto max_cost = std::max_element(cost.begin(), cost.end());
    auto min_time = std::min_element(time.begin(), time.end());
    auto r0 = max_cost * G.number_of_edges() / min_time;
    return max_parametric(G, r0, calc_weight, calc_ratio);
}

// if (__name__ == "__main__") {
//     #include <xnetwork.hpp> // as xn
//     // from neg_cycle import *
//     // from xnetwork.utils import generate_unique_node

//     G = create_test_case1();
//     G[1][2]["cost"] = 5
//     r, c, dist = min_cycle_ratio(G);
//     CHECK(c != None
//     print(r);
//     print(c);
//     print(dist.items());

//     G = xn::cycle_graph(5, create_using=xn::DiGraph());
//     G[1][2]["cost"] = -6.
//     newnode = generate_unique_node();
//     G.add_edges_// from([(newnode, n) for (auto n : G]);
//     r, c, dist = min_cycle_ratio(G);
//     CHECK(c != None
//     print(r);
//     print(c);
//     print(dist.items());

#endif