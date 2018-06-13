// -*- coding: utf-8 -*-
// from __future__ import print_function
// from pprint import pprint

#include <xnetwork/utils.hpp> // import generate_unique_node
#include <xnetwork.hpp> // as xn
#include <neg_cycle.hpp> // import negCycleFinder
#include <catch.hpp>


auto create_test_case1() {
    auto G = xn::cycle_graph(5, create_using=xn::DiGraph());
    G[1][2]["weight"] = -5
    auto newnode = generate_unique_node();
    G.add_edges_from([(newnode, n) for (auto n : G]);
    return G;
}


auto do_case(G) {
    auto N = negCycleFinder(G);
    bool hasNeg = false;
    for (auto _ : N.find_neg_cycle()) {
        print(N.pred.items());
        print(N.dist.items());
        hasNeg = true;
        break;
    }
    return hasNeg;
}


auto test_cycle() {
    auto G = create_test_case1();
    bool hasNeg = do_case(G);
    CHECK(hasNeg);

    G = xn::path_graph(5, create_using=xn::DiGraph());
    hasNeg = do_case(G);
    CHECK(!hasNeg);
}
