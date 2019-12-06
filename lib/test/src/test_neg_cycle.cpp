// -*- coding: utf-8 -*-
#include <array>
#include <catch2/catch.hpp>
#include <netoptim/neg_cycle.hpp> // import negCycleFinder
#include <vector>
#include <xnetwork/classes/digraphs.hpp>

/*!
 * @brief Create a test case1 object
 *
 * @return auto
 */
static auto create_test_case1()
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
    auto edges =
        std::array {Edge(A, B), Edge(B, C), Edge(C, D), Edge(D, E), Edge(E, A)};
    constexpr auto weights = std::array {-5, 1, 1, 1, 1};
    auto g = xn::DiGraphS(py::range<int>(num_nodes));
    g.add_edges_from(edges, weights);
    return g;
}

/*!
 * @brief Create a test case2 object
 *
 * @return auto
 */
static auto create_test_case2()
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
    auto edges =
        std::array {Edge(A, B), Edge(B, C), Edge(C, D), Edge(D, E), Edge(E, A)};
    auto weights = std::array {2, 1, 1, 1, 1};
    auto g = xn::DiGraphS(py::range<int>(num_nodes));
    g.add_edges_from(edges, weights);
    return g;
}

/*!
 * @brief Create a test case timing object
 *
 * @return auto
 */
static auto create_test_case_timing()
{
    using Edge = std::pair<int, int>;
    constexpr auto num_nodes = 3;
    enum nodes
    {
        A,
        B,
        C
    };
    const auto edges = std::array {Edge(A, B), Edge(B, A), Edge(B, C),
        Edge(C, B), Edge(B, C), Edge(C, B), Edge(C, A), Edge(A, C)};
    const auto weights = std::array {7, 0, 3, 1, 6, 4, 2, 5};
    auto g = xn::DiGraphS(py::range<int>(num_nodes));
    g.add_edges_from(edges, weights);
    return g;
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
bool do_case(const Graph& G)
{
    auto get_weight = [&](const auto& e) -> int {
        const auto [u, v] = G.end_points(e);
        return G[u][v];
    };

    auto dist = std::vector(G.number_of_nodes(), 0);
    auto N = negCycleFinder(G);
    auto cycle = N.find_neg_cycle(dist, get_weight);
    return !cycle.empty();
}

/*!
 * @brief
 *
 */
TEST_CASE("Test Negative Cycle", "[test_neg_cycle]")
{
    auto G = create_test_case1();
    auto hasNeg = do_case(G);
    CHECK(hasNeg);
}

/*!
 * @brief
 *
 */
TEST_CASE("Test No Negative Cycle", "[test_neg_cycle]")
{
    auto G = create_test_case2();
    auto hasNeg = do_case(G);
    CHECK(!hasNeg);
}

/*!
 * @brief
 *
 */
TEST_CASE("Test Timing Graph", "[test_neg_cycle]")
{
    auto G = create_test_case_timing();
    auto hasNeg = do_case(G);
    CHECK(!hasNeg);
}
