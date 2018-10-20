// -*- coding: utf-8 -*-
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>
#include <catch.hpp>
#include <ellcpp/oracles/fraction.hpp>        // import Fraction
#include <ellcpp/oracles/min_cycle_ratio.hpp> // import min_cycle_ratio, set_default
#include <py2cpp/nx2bgl.hpp>
#include <utility> // for std::pair
#include <vector>
// from fractions import Fraction

namespace boost {

enum edge_id_tag_t { id_tag }; // a unique #
BOOST_INSTALL_PROPERTY(edge, id_tag);

} // namespace boost

using graph_t = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS, boost::no_property,
    boost::property<boost::edge_id_tag_t, std::size_t>>;
using Vertex = typename boost::graph_traits<graph_t>::vertex_descriptor;
using edge_t = typename boost::graph_traits<graph_t>::edge_iterator;

static auto create_test_case1() {
    using Edge = std::pair<int, int>;
    const int num_nodes = 5;
    enum nodes { A, B, C, D, E };
    // char name[] = "ABCDE";
    static Edge edge_array[] = {Edge(A, B), Edge(B, C), Edge(C, D), Edge(D, E),
                                Edge(E, A)};
    std::size_t indices[] = {0, 1, 2, 3, 4};
    int num_arcs = sizeof(edge_array) / sizeof(Edge);
    static graph_t g(edge_array, edge_array + num_arcs, indices, num_nodes);
    return xn::grAdaptor<graph_t>(g);
}

static auto create_test_case_timing() {
    using Edge = std::pair<int, int>;
    const int num_nodes = 3;
    enum nodes { A, B, C };
    // char name[] = "ABCDE";
    static Edge edge_array[] = {Edge(A, B), Edge(B, A), Edge(B, C), Edge(C, B),
                                Edge(B, C), Edge(C, B), Edge(C, A), Edge(A, C)};
    std::size_t indices[] = {0, 1, 2, 3, 4, 5, 6, 7};
    // int weights[] = {7, 0, 3, 1, 2, 5};
    int num_arcs = sizeof(edge_array) / sizeof(Edge);
    static graph_t g(edge_array, edge_array + num_arcs, indices, num_nodes);
    return xn::grAdaptor<graph_t>(g);
}


TEST_CASE("Test Cycle Ratio", "[test_cycle_ratio]") {
    using EdgeIndexMap =
        typename boost::property_map<graph_t, boost::edge_id_tag_t>::type;
    using IterMap =
        boost::iterator_property_map<fun::Fraction<int> *, EdgeIndexMap,
                                     fun::Fraction<int>, fun::Fraction<int> &>;

    xn::grAdaptor<graph_t> G = create_test_case1();
    fun::Fraction<int> cost[] = {{5, 1}, {1, 1}, {1, 1}, {1, 1}, {1, 1}};
    fun::Fraction<int> time[] = {{1, 1}, {1, 1}, {1, 1}, {1, 1}, {1, 1}};
    EdgeIndexMap edge_id = boost::get(boost::id_tag, G);
    IterMap cost_pa(cost, edge_id), time_pa(time, edge_id);

    auto get_cost = [&](const xn::grAdaptor<graph_t>& G, const auto& e) -> fun::Fraction<int> {
        return boost::get(cost_pa, e);
    };
    auto get_time = [&](const xn::grAdaptor<graph_t>& G, const auto& e) -> fun::Fraction<int> {
        return boost::get(time_pa, e);
    };

    auto [r, c] = min_cycle_ratio(G, get_cost, get_time, fun::Fraction<int>{});
    CHECK(!c.empty());
    CHECK(c.size() == 5);
    CHECK(r == fun::Fraction<int>(9, 5));
    // print(r);
    // print(c);
    // print(dist.items());
}


TEST_CASE("Test Cycle Ratio of Timing Graph", "[test_cycle_ratio]") {
    using EdgeIndexMap =
        typename boost::property_map<graph_t, boost::edge_id_tag_t>::type;
    using IterMap =
        boost::iterator_property_map<fun::Fraction<int> *, EdgeIndexMap,
                                     fun::Fraction<int>, fun::Fraction<int> &>;

    xn::grAdaptor<graph_t> G = create_test_case_timing();
    fun::Fraction<int> cost[] = {{7, 1}, {0, 1}, {3, 1}, {1, 1}, {6, 1}, {4, 1}, {2, 1}, {5, 1}};
    fun::Fraction<int> time[] = {{1, 1}, {1, 1}, {1, 1}, {1, 1}, {1, 1}, {1, 1}, {1, 1}, {1, 1}};
    EdgeIndexMap edge_id = boost::get(boost::id_tag, G);
    IterMap cost_pa(cost, edge_id), time_pa(time, edge_id);

    auto get_cost = [&](const xn::grAdaptor<graph_t>& G, const auto& e) -> fun::Fraction<int> {
        return boost::get(cost_pa, e);
    };
    auto get_time = [&](const xn::grAdaptor<graph_t>& G, const auto& e) -> fun::Fraction<int> {
        return boost::get(time_pa, e);
    };

    auto [r, c] = min_cycle_ratio(G, get_cost, get_time, fun::Fraction<int>{});
    CHECK(!c.empty());
    CHECK(r == fun::Fraction<int>(2, 1));
    CHECK(c.size() == 2);
    // print(r);
    // print(c);
    // print(dist.items());
}
