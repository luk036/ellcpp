// -*- coding: utf-8 -*-

#include <catch.hpp>
#include <utility>                   // for std::pair
#include <py2cpp/nx2bgl.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <ellcpp/oracles/min_cycle_ratio.hpp> // import min_cycle_ratio, set_default
// from fractions import Fraction

typedef boost::adjacency_list < boost::listS, boost::vecS, boost::directedS, boost::no_property, boost::property < boost::edge_weight_t, int > > graph_t;
typedef boost::graph_traits< graph_t >::vertex_descriptor Vertex;
//typedef boost::graph_traits< graph_t >::edge_iterator edge_t;


static auto create_test_case1() {
    using Edge = std::pair<int, int>;
    const int num_nodes = 5;
    enum nodes { A, B, C, D, E };
    char name[] = "ABCDE";
    static Edge edge_array[] = { Edge(A, B), Edge(B, C), Edge(C, D), Edge(D, E), Edge(E, A) };
    int weights[] = { -5, 1, 1, 1, 1 };
    int num_arcs = sizeof(edge_array) / sizeof(Edge);
    static graph_t g(edge_array, edge_array + num_arcs, weights, num_nodes);
    return xn::grAdaptor<graph_t>(g);
}


TEST_CASE("Test Cycle Ratio", "[test_cycle_ratio]") {
    xn::grAdaptor<graph_t> G = create_test_case1();
    using edge_t = decltype(*(G.edges().begin()));
    py::dict<edge_t*, double> cost;
    py::dict<edge_t*, double> time;

    set_default(G, cost, 1.);
    cost[&G.edges().begin()] = 5.;
    auto [r, c, dist] = min_cycle_ratio(G, cost, time);
    CHECK(!c.empty());
    CHECK(r == 9./5.);
    // print(r);
    // print(c);
    // print(dist.items());
}
