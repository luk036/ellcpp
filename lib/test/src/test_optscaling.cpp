// -*- coding: utf-8 -*-
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>
#include <catch.hpp>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/optscaling_oracle.hpp> // import optscaling
#include <py2cpp/nx2bgl.hpp>
#include <utility> // for std::pair
// #include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

namespace boost
{

enum edge_id_tag_t
{
    id_tag
}; // a unique #
BOOST_INSTALL_PROPERTY(edge, id_tag);

} // namespace boost

using Arr = xt::xarray<double, xt::layout_type::row_major>;

using graph_t =
    boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property,
                          boost::property<boost::edge_id_tag_t, std::size_t>>;
using Vertex = boost::graph_traits<graph_t>::vertex_descriptor;
using edge_t = boost::graph_traits<graph_t>::edge_iterator;

static xn::grAdaptor<graph_t> create_test_case1()
{
    using Edge          = std::pair<int, int>;
    const auto num_nodes = 5;
    enum nodes
    {
        A,
        B,
        C,
        D,
        E
    };
    // char name[] = "ABCDE";
    static Edge    edge_array[] = {Edge(A, B), Edge(B, C), Edge(C, D), Edge(D, E), Edge(E, A)};
    std::size_t    indices[]    = {0, 1, 2, 3, 4};
    int            num_arcs     = sizeof(edge_array) / sizeof(Edge);
    static graph_t g(edge_array, edge_array + num_arcs, indices, num_nodes);
    return xn::grAdaptor<graph_t>(g);
}

TEST_CASE("Test Optimal Scaling", "[test_optscaling]")
{
    using EdgeIndexMap = typename boost::property_map<graph_t, boost::edge_id_tag_t>::type;
    using IterMap      = boost::iterator_property_map<double*, EdgeIndexMap, double, double&>;

    auto G = create_test_case1();

    double    elem[]       = {1.2, 2.3, 3.4, -4.5, 5.6};
    const auto num_of_nodes = sizeof(elem) / sizeof(double);
    double    cost[num_of_nodes];
    for (auto i = 0; i < num_of_nodes; ++i)
    {
        cost[i] = std::log(std::abs(elem[i]));
    }
    auto edge_id = boost::get(boost::id_tag, G);
    auto cost_pa = IterMap{cost, edge_id};

    auto get_cost = [&](const xn::grAdaptor<graph_t>& G, const auto& e) -> double {
        return boost::get(cost_pa, e);
    };

    auto [min, max] = std::minmax_element(cost, cost + num_of_nodes);
    // auto cmin = *std::min_element(cost, cost + num_of_nodes);
    auto x0       = Arr{*max, *min};
    auto t        = *max - *min;
    auto E        = ell{1.5 * t, x0};
    auto P        = optscaling_oracle{G, get_cost, double(0.)};
    auto ell_info = cutting_plane_dc(P, E, std::numeric_limits<double>::max());

    CHECK(ell_info.feasible);
    CHECK(ell_info.num_iters <= 27);
}
