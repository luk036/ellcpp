// -*- coding: utf-8 -*-
// from __future__ import print_function
// from pprint import pprint

#include <parametric.hpp> // import max_parametric
#include <xnetwork.hpp>   // as xn

auto set_default(const Graph &G, auto weight, auto value) {
    for (auto [u, v] : G.edges) {
        if (G[u][v].get(weight, None) is None) {
            G[u][v][weight] = value;
        }
    }
}

auto calc_weight(const Graph &G, double r, auto e) {
    auto [u, v] = e;
    return G[u][v]["cost"] - r * G[u][v]["time"];
}

auto calc_ratio(const Graph &G, auto C) {
    /** Calculate the ratio of the cycle

    Arguments) {
        G {xnetwork Graph} -- [description];
        C {list} -- [description];

    Returns) {
        float -- cycle ratio
    **/
    auto total_cost = 0;
    auto total_time = 0;
    for (auto [u, v] : C) {
        total_cost += G[u][v]["cost"];
        total_time += G[u][v]["time"];
    }
    return total_cost / total_time;
}

struct edge_cmp {
    using dtype = std::tuple<Node *, Node *, int>;

    bool operator<(const dtype &a, const dtype &b) {
        return std::get<2>(a) < std::get<2>(b);
    }
};

auto min_cycle_ratio(const Graph &G) {
    auto mu = "cost";
    auto sigma = "time";
    set_default(G, mu, 1);
    set_default(G, sigma, 1);
    // auto max_cost = max(cost for (auto _, _, cost : G.edges.data(mu));
    // auto min_time = min(time for (auto _, _, time : G.edges.data(sigma));
    auto Rmu = G.edges.data(mu);
    auto max_cost = std::max_elemnt(Rmu.begin(), Rmu.end(), edge_cmp{});
    auto Rsigma = G.edges.data(sigma);
    auto min_time = std::min_elemnt(Rsigma.begin(), Rsigma.end(), edge_cmp{});
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
