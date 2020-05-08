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
    const auto get_weight = [&](const auto& e) -> int {
        const auto [u, v] = G.end_points(e);
        return G[u].at(v);
    };

    auto dist = std::vector(G.number_of_nodes(), 0);
    auto N = negCycleFinder(G);
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
    const auto edges = std::array {
        Edge {A, B}, Edge {B, C}, Edge {C, D}, Edge {D, E}, Edge {E, A}};
    constexpr auto weights = std::array {-5, 1, 1, 1, 1};
    auto G = xn::DiGraphS {py::range<int>(num_nodes)};
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
    auto edges = std::array {
        Edge {A, B}, Edge {B, C}, Edge {C, D}, Edge {D, E}, Edge {E, A}};
    auto weights = std::array {2, 1, 1, 1, 1};
    auto G = xn::DiGraphS {py::range<int>(num_nodes)};
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
    const auto edges = std::array {Edge {A, B}, Edge {B, A}, Edge {B, C},
        Edge {C, B}, Edge {B, C}, Edge {C, B}, Edge {C, A}, Edge {A, C}};
    const auto weights = std::array {7, 0, 3, 1, 6, 4, 2, 5};
    auto G = xn::DiGraphS {py::range<int>(num_nodes)};
    G.add_edges_from(edges, weights);
    const auto hasNeg = do_case(G);
    CHECK(!hasNeg);
}
