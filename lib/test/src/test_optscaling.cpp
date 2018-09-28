// -*- coding: utf-8 -*-
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>
#include <catch.hpp>
#include <py2cpp/nx2bgl.hpp>
#include <utility> // for std::pair
#include <vector>
// from fractions import Fraction
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/optscaling_oracle.hpp> // import optscaling
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <algorithm>

namespace boost {

enum edge_id_tag_t { id_tag }; // a unique #
BOOST_INSTALL_PROPERTY(edge, id_tag);

} // namespace boost

using Arr = xt::xarray<double>;

typedef boost::adjacency_list<
    boost::vecS, boost::vecS, boost::directedS, boost::no_property,
    boost::property<boost::edge_id_tag_t, std::size_t>>
    graph_t;
typedef boost::graph_traits<graph_t>::vertex_descriptor Vertex;
typedef boost::graph_traits<graph_t>::edge_iterator edge_t;

static auto create_test_case1() {
    using Edge = std::pair<int, int>;
    const int num_nodes = 5;
    enum nodes { A, B, C, D, E };
    char name[] = "ABCDE";
    static Edge edge_array[] = {Edge(A, B), Edge(B, C), Edge(C, D), Edge(D, E),
                                Edge(E, A)};
    std::size_t indices[] = {0, 1, 2, 3, 4};
    int num_arcs = sizeof(edge_array) / sizeof(Edge);
    static graph_t g(edge_array, edge_array + num_arcs, indices, num_nodes);
    return xn::grAdaptor<graph_t>(g);
}

TEST_CASE("Test Optimal Scaling", "[test_optscaling]") {
    using EdgeIndexMap =
        typename boost::property_map<graph_t, boost::edge_id_tag_t>::type;
    using IterMap =
        boost::iterator_property_map<double *, EdgeIndexMap, double, double &>;

    xn::grAdaptor<graph_t> G = create_test_case1();
    const int N = 5;
    double elem[N] = {1.2, 2.3, 3.4, -4.5, 5.6};
    double cost[N];
    for (auto i = 0; i < N; ++i) {
        cost[i] = std::log(std::abs(elem[i]));
    }
    EdgeIndexMap edge_id = boost::get(boost::id_tag, G);
    IterMap cost_pa(cost, edge_id);

    auto cmax = *std::max_element(cost, cost+N);
    auto cmin = *std::min_element(cost, cost+N);

    auto x0 = Arr{cmax, cmin};
    auto t = cmax - cmin;
    auto E = ell(1.5*t, x0);
    auto P = optscaling_oracle(G, cost_pa, double(0.0));
    auto [xb, fb, niter, feasible, status] = 
        cutting_plane_dc(P, E, 1.001*t);
    CHECK(feasible);
    CHECK(niter <= 27);
}


