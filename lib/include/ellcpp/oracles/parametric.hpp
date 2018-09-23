// -*- coding: utf-8 -*-
// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_PARAMETRIC_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_PARAMETRIC_HPP 1

#include "neg_cycle.hpp" // import negCycleFinder
#include <tuple>
#include <vector>

/** 
 * maximum parametric problem:
 *    max  r
 *    s.t. dist[v] - dist[v] <= d(u,v,r);
 *         for (auto all (u, v] : G
 *
 *   Arguments:
 *       G {[type]} -- [description];
 *       r {float} -- parameter to be maximized, initially a large number (infeasible) 
 *       d {[type]} -- monotone decreasing function w.r.t. r 
 *       zero_cancel {[type]} -- [description];
 *
 *   Returns:
 *       r_opt -- optimal value
 *       C_opt -- Most critial cycle
 *       dist -- optimal sol"n
 */
template <typename Graph, typename Fn1, typename Fn2>
auto max_parametric(const Graph &G, double r, Fn1 &d, Fn2 &zero_cancel) {
    auto get_weight = [](const Graph &G, auto &e) -> int { // int???
        return d(G, r, e); 
    };

    auto S = negCycleFinder(G, get_weight, double());
    std::vector<std::tuple<Node *, Node *>> C_opt{};
    auto r_opt = r;

    while (true) {
        auto C = S.neg_cycle_relax();
        auto r_min = zero_cancel(G, C);

        if (r_min >= r_opt) {
            break;
        }
        C_opt = C;
        r_opt = r_min;
        // update ???
        for (auto e : C_opt) {
            auto u = G.source(e);
            auto v = G.target(e);
            S.dist[u] = S.dist[v] - get_weight(G, e);
        }
    }

    return std::tuple{r_opt, C_opt, S._dist};
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