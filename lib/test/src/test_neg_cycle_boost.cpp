// -*- coding: utf-8 -*-
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <catch2/catch.hpp>
#include <netoptim/min_cycle_ratio.hpp>
#include <netoptim/neg_cycle.hpp> // import negCycleFinder
#include <py2cpp/nx2bgl.hpp>
#include <utility> // for std::pair

using graph_t = boost::adjacency_list<boost::listS, boost::vecS,
    boost::directedS, boost::no_property,
    boost::property<boost::edge_weight_t, int,
        boost::property<boost::edge_index_t, int>>>;
using Vertex = boost::graph_traits<graph_t>::vertex_descriptor;
using Edge_it = boost::graph_traits<graph_t>::edge_iterator;

static xn::grAdaptor<graph_t> create_test_case1()
{
    using Edge = std::pair<int, int>;
    const auto num_nodes = 5;
    enum nodes
    {
        A,
        B,
        C,
        D,
        E
    };
    static Edge edge_array[] = {
        Edge {A, B}, Edge {B, C}, Edge {C, D}, Edge {D, E}, Edge {E, A}};
    int weights[] = {-5, 1, 1, 1, 1};
    int num_arcs = sizeof(edge_array) / sizeof(Edge);
    static auto g =
        graph_t(edge_array, edge_array + num_arcs, weights, num_nodes);
    return xn::grAdaptor<graph_t> {g};
}

static xn::grAdaptor<graph_t> create_test_case2()
{
    using Edge = std::pair<int, int>;
    const auto num_nodes = 5;
    enum nodes
    {
        A,
        B,
        C,
        D,
        E
    };
    static Edge edge_array[] = {
        Edge {A, B}, Edge {B, C}, Edge {C, D}, Edge {D, E}, Edge {E, A}};
    int weights[] = {2, 1, 1, 1, 1};
    int num_arcs = sizeof(edge_array) / sizeof(Edge);
    static auto g =
        graph_t(edge_array, edge_array + num_arcs, weights, num_nodes);
    return xn::grAdaptor<graph_t> {g};
}

static auto create_test_case_timing() -> xn::grAdaptor<graph_t>
{
    using Edge = std::pair<int, int>;
    constexpr auto num_nodes = 3;
    enum nodes
    {
        A,
        B,
        C
    };
    static Edge edge_array[] = {Edge {A, B}, Edge {B, A}, Edge {B, C},
        Edge {C, B}, Edge {B, C}, Edge {C, B}, Edge {C, A}, Edge {A, C}};
    int weights[] = {7, 0, 3, 1, 6, 4, 2, 5};
    constexpr int num_arcs = sizeof(edge_array) / sizeof(Edge);
    static auto g =
        graph_t(edge_array, edge_array + num_arcs, weights, num_nodes);
    return xn::grAdaptor<graph_t> {g};
}

auto do_case(const xn::grAdaptor<graph_t>& G) -> bool
{
    using edge_t = decltype(*(std::begin(G.edges())));

    const auto get_weight = [&](const edge_t& e) -> int {
        const auto& weightmap = boost::get(boost::edge_weight, G);
        return weightmap[e];
    };

    auto dist = std::vector(G.number_of_nodes(), 0);
    auto N = negCycleFinder {G};
    const auto cycle = N.find_neg_cycle(dist, get_weight);
    return !cycle.empty();
}

TEST_CASE("Test Negative Cycle (boost)", "[test_neg_cycle_boost]")
{
    const auto G = create_test_case1();
    const auto hasNeg = do_case(G);
    CHECK(hasNeg);
}

TEST_CASE("Test No Negative Cycle (boost)", "[test_neg_cycle_boost]")
{
    const auto G = create_test_case2();
    const auto hasNeg = do_case(G);
    CHECK(!hasNeg);
}

TEST_CASE("Test Timing Graph (boost)", "[test_neg_cycle_boost]")
{
    const auto G = create_test_case_timing();
    const auto hasNeg = do_case(G);
    CHECK(!hasNeg);
}
