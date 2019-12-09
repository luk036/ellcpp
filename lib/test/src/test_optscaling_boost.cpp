// -*- coding: utf-8 -*-
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>
#include <catch2/catch.hpp>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
// #include <ellcpp/ell1d.hpp>
#include <ellcpp/oracles/optscaling_oracle.hpp> // import optscaling
#include <py2cpp/nx2bgl.hpp>
#include <utility> // for std::pair
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
    boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
        boost::no_property, boost::property<boost::edge_id_tag_t, size_t>>;
using Vertex = typename boost::graph_traits<graph_t>::vertex_descriptor;
using edge_t = typename boost::graph_traits<graph_t>::edge_iterator;

/*!
 * @brief Create a test case1 object
 *
 * @return xn::grAdaptor<graph_t>
 */
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
    // char name[] = "ABCDE";
    Edge edge_array[] = {
        Edge(A, B), Edge(B, C), Edge(C, D), Edge(D, E), Edge(E, A)};
    size_t indices[] = {0, 1, 2, 3, 4};
    auto num_arcs = sizeof(edge_array) / sizeof(Edge);
    static auto g =
        graph_t(edge_array, edge_array + num_arcs, indices, num_nodes);
    return xn::grAdaptor<graph_t>(g);
}

TEST_CASE(
    "Test Optimal Scaling (two varaibles, boost)", "[test_optscaling_boost]")
{
    using EdgeIndexMap =
        typename boost::property_map<graph_t, boost::edge_id_tag_t>::type;
    using IterMap =
        boost::iterator_property_map<double*, EdgeIndexMap, double, double&>;

    auto G = create_test_case1();

    double elem[] = {1.2, 2.3, 3.4, -4.5, 5.6};
    const auto num_of_nodes = sizeof(elem) / sizeof(double);

    double cost[num_of_nodes];
    for (size_t i = 0U; i != num_of_nodes; ++i)
    {
        cost[i] = std::log(std::abs(elem[i]));
    }
    auto edge_id = boost::get(boost::id_tag, G);
    auto cost_pa = IterMap {cost, edge_id};

    auto get_cost = [&](const auto& e) -> double {
        return boost::get(cost_pa, e);
    };

    const auto [cmin, cmax] =
        std::minmax_element(std::begin(cost), std::end(cost));
    // auto cmin = *std::min_element(cost, cost + num_of_nodes);
    const auto x0 = Arr {*cmax, *cmin};
    auto t1 = *cmax - *cmin;
    auto E = ell {1.5 * t1, x0};
    auto dist = std::vector(G.number_of_nodes(), 0.);

    auto P = optscaling_oracle {G, dist, get_cost};
    auto t = std::numeric_limits<double>::max();
    const auto [x, ell_info] = cutting_plane_dc(P, E, t);

    CHECK(x[1] >= x[0]);
    CHECK(ell_info.feasible);
    CHECK(ell_info.num_iters <= 27);
}
