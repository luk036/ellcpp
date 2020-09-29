// -*- coding: utf-8 -*-
#include <array>
#include <doctest.h>
#include <netoptim/neg_cycle.hpp> // import negCycleFinder
#include <vector>
#include <xnetwork/classes/digraphs.hpp>


/*!
 * @brief
 *
 * @tparam Graph
 * @param G
 * @return true
 * @return false
 */
template <typename Graph>
bool do_case(const Graph& G)
{
    const auto get_weight = [&](const auto& edge) -> int {
        const auto e = G.end_points(edge);
        return G[e.first].at(e.second);
    };

    auto dist = std::vector<int>(G.number_of_nodes(), 0);
    auto N = negCycleFinder<Graph>(G);
    const auto cycle = N.find_neg_cycle(dist, get_weight);
    return !cycle.empty();
}

/*!
 * @brief
 *
 */
TEST_CASE("Test Negative Cycle")
{
    using Edge = std::pair<int, int>;
    constexpr int num_nodes = 5;
    enum nodes
    {
        A,
        B,
        C,
        D,
        E
    };
    const auto edges = std::array<Edge, 5> {
        Edge {A, B}, Edge {B, C}, Edge {C, D}, Edge {D, E}, Edge {E, A}};
    constexpr auto weights = std::array<int, 5> {-5, 1, 1, 1, 1};
    auto G = xn::SimpleDiGraphS {py::range<int>(num_nodes)};
    G.add_edges_from(edges, weights);
    const auto hasNeg = do_case(G);
    CHECK(hasNeg);
}

/*!
 * @brief
 *
 */
TEST_CASE("Test No Negative Cycle")
{
    using Edge = std::pair<int, int>;
    constexpr auto num_nodes = 5;
    enum nodes
    {
        A,
        B,
        C,
        D,
        E
    };
    auto edges = std::array<Edge, num_nodes> {
        Edge {A, B}, Edge {B, C}, Edge {C, D}, Edge {D, E}, Edge {E, A}};
    auto weights = std::array<int, num_nodes> {2, 1, 1, 1, 1};
    auto G = xn::SimpleDiGraphS {py::range<int>(num_nodes)};
    G.add_edges_from(edges, weights);
    const auto hasNeg = do_case(G);
    CHECK(!hasNeg);
}

/*!
 * @brief
 *
 */
TEST_CASE("Test Timing Graph")
{
    using Edge = std::pair<int, int>;
    constexpr auto num_nodes = 3;
    enum nodes
    {
        A,
        B,
        C
    };
    const auto edges =
        std::array<Edge, 8> {Edge {A, B}, Edge {B, A}, Edge {B, C}, Edge {C, B},
            Edge {B, C}, Edge {C, B}, Edge {C, A}, Edge {A, C}};
    const auto weights = std::array<int, 8> {7, 0, 3, 1, 6, 4, 2, 5};
    auto G = xn::SimpleDiGraphS {py::range<int>(num_nodes)};
    G.add_edges_from(edges, weights);
    const auto hasNeg = do_case(G);
    CHECK(!hasNeg);
}

/*!
 * @brief
 *
 */
TEST_CASE("Test Timing Graph (2)")
{
    using Edge = std::pair<int, int>;
    constexpr auto num_nodes = 3;
    enum nodes
    {
        A,
        B,
        C
    };
    const auto edges =
        std::array<Edge, 8> {Edge {A, B}, Edge {B, A}, Edge {B, C}, Edge {C, B},
            Edge {B, C}, Edge {C, B}, Edge {C, A}, Edge {A, C}};
    const auto weights = std::array<int, 8> {3, -4, -1, -1, 2, 0, -2, 1};
    auto G = xn::SimpleDiGraphS {py::range<int>(num_nodes)};
    G.add_edges_from(edges, weights);
    const auto hasNeg = do_case(G);
    CHECK(hasNeg);
}


/*!
 * @brief
 *
 * @tparam Graph
 * @param G
 * @return true
 * @return false
 */
template <typename Graph>
bool do_case_float(const Graph& G)
{
    const auto get_weight = [&](const auto& edge) -> double {
        const auto e = G.end_points(edge);
        return G[e.first].at(e.second);
    };

    auto dist = std::vector<double>(G.number_of_nodes(), 0.0);
    auto N = negCycleFinder<Graph>(G);
    const auto cycle = N.find_neg_cycle(dist, get_weight);
    return !cycle.empty();
}

/*!
 * @brief
 *
 */
TEST_CASE("Test Negative Cycle")
{
    using Edge = std::pair<int, int>;
    constexpr int num_nodes = 5;
    enum nodes
    {
        A,
        B,
        C,
        D,
        E
    };
    const auto edges = std::array<Edge, 5> {
        Edge {A, B}, Edge {B, C}, Edge {C, D}, Edge {D, E}, Edge {E, A}};
    constexpr auto weights = std::array<int, 5> {-5, 1, 1, 1, 1};
    auto G = xn::SimpleDiGraphS {py::range<int>(num_nodes)};
    G.add_edges_from(edges, weights);
    const auto hasNeg = do_case(G);
    CHECK(hasNeg);
}
