// -*- coding: utf-8 -*-
#include <catch2/catch.hpp>
#include <utility> // for std::pair

// auto get_weight(const graph_t &G, const Edge_it &e) {
//     auto weightmap = boost::get(boost::edge_weight_t(), G);
//     return weightmap[*e];
//     //auto u = boost::source(e, G);
//     //auto v = boost::target(e, G);
//     //return G[u][v].get("weight", 1);
// }

// static graph_t create_test_case1() {
//     using Edge = std::pair<int, int>;
//     const auto num_nodes = 5;
//     enum nodes { A, B, C, D, E };
//     // char name[] = "ABCDE";
//     static Edge edge_array[] = {Edge(A, B), Edge(B, C), Edge(C, D), Edge(D,
//     E),
//                                 Edge(E, A)};
//     int weights[] = {-5, 1, 1, 1, 1};
//     int num_arcs = sizeof(edge_array) / sizeof(Edge);
//     graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);
//     return std::move(g);
// }

// static graph_t create_test_case2() {
//     using Edge = std::pair<int, int>;
//     const auto num_nodes = 5;
//     enum nodes { A, B, C, D, E };
//     // char name[] = "ABCDE";
//     static Edge edge_array[] = {Edge(A, B), Edge(B, C), Edge(C, D), Edge(D,
//     E),
//                                 Edge(E, A)};
//     int weights[] = {2, 1, 1, 1, 1};
//     int num_arcs = sizeof(edge_array) / sizeof(Edge);
//     graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);
//     return std::move(g);
// }

// static graph_t create_test_case_timing() {
//     using Edge = std::pair<int, int>;
//     const auto num_nodes = 3;
//     enum nodes { A, B, C };
//     // char name[] = "ABCDE";
//     static Edge edge_array[] = {Edge(A, B), Edge(B, A), Edge(B, C), Edge(C,
//     B),
//                                 Edge(B, C), Edge(C, B), Edge(C, A), Edge(A,
//                                 C)};
//     int weights[] = {7, 0, 3, 1, 6, 4, 2, 5};
//     int num_arcs = sizeof(edge_array) / sizeof(Edge);
//     graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);
//     return std::move(g);
// }

// TEST_CASE("Test XNetwork", "[test_xnetwork]") {
//     auto G = create_test_case1();
//     // CHECK(hasNeg);
// }

// TEST_CASE("Test XNetwork 2", "[test_xnetwork]") {
//     auto G = create_test_case2();
//     // CHECK(!hasNeg);
// }

// TEST_CASE("Test Timing Graph", "[test_xnetwork]") {
//     auto G = create_test_case_timing();
//     // CHECK(!hasNeg);
// }
