// -*- coding: utf-8 -*-
#pragma once

/*!
Negative cycle detection for (auto weighed graphs.
**/
#include <cassert>
#include <py2cpp/py2cpp.hpp>
#include <vector>

/*!
 * @brief negative cycle
 *
 * @tparam Graph
 *
 * Note: Bellman-Ford's shortest-path algorithm (BF) is NOT the best way to
 *       detect negative cycles, because
 *
 *  1. BF needs a source node.
 *  2. BF detect whether there is a negative cycle at the fianl stage.
 *  3. BF restarts the solution (dist[u]) every time.
 */
template <typename Graph> //
class negCycleFinder
{
    using node_t = typename Graph::node_t;
    using edge_t = typename Graph::edge_t;

    py::dict<node_t, node_t> _pred {};
    py::dict<node_t, edge_t> _edge {};

  private:
    const Graph& _G; // const???

  public:
    /*!
     * @brief Construct a new neg Cycle Finder object
     *
     * @param G
     */
    explicit negCycleFinder(const Graph& G)
        : _G {G}
    {
    }

  public:
    /*!
     * @brief find negative cycle
     *
     * @tparam Container
     * @tparam WeightFn
     * @param dist
     * @param get_weight
     * @return std::vector<edge_t>
     */
    template <typename Container, typename WeightFn>
    auto find_neg_cycle(Container& dist, const WeightFn& get_weight)
        -> std::vector<edge_t>
    {
        this->_pred.clear();
        this->_edge.clear();

        while (true)
        {
            auto changed = this->__relax(dist, get_weight);
            if (!changed)
            {
                break;
            }
            auto v = this->__find_cycle();
            if (v != this->_G.null_vertex())
            {
                assert(this->__is_negative(v, dist, get_weight));
                return this->__cycle_list(v);
            }
        }
        return std::vector<edge_t> {}; // ???
    }

  private:
    /*!
     * @brief Find a cycle on policy graph
     *
     * @return node_t a start node of the cycle
     */
    auto __find_cycle() -> node_t
    {
        auto visited = py::dict<node_t, node_t> {};

        for (auto v : this->_G)
        {
            if (visited.contains(v))
            {
                continue;
            }
            auto u = v;
            while (true)
            {
                visited[u] = v;
                if (!this->_pred.contains(u))
                {
                    break;
                }
                u = this->_pred[u];
                if (visited.contains(u))
                {
                    if (visited[u] == v)
                    {
                        // if (this->__is_negative(u)) {
                        // should be "yield u";
                        return u;
                        // }
                    }
                    break;
                }
            }
        }

        return this->_G.null_vertex();
    }

    /*!
     * @brief Perform one relaxation
     *
     * @tparam Container
     * @tparam WeightFn
     * @param dist
     * @param get_weight
     * @return true
     * @return false
     */
    template <typename Container, typename WeightFn>
    auto __relax(Container& dist, const WeightFn& get_weight) -> bool
    {
        auto changed = false;
        for (auto e : this->_G.edges())
        {
            auto [u, v] = this->_G.end_points(e);

            auto wt = get_weight(this->_G, e);
            auto d = dist[u] + wt;
            if (dist[v] > d)
            {
                dist[v] = d;
                this->_pred[v] = u;
                this->_edge[v] = e; // ???
                changed = true;
            }
        }
        return changed;
    }

    /*!
     * @brief generate a cycle list
     *
     * @param handle
     * @return std::vector<edge_t>
     */
    auto __cycle_list(node_t handle) -> std::vector<edge_t>
    {
        auto v = handle;
        auto cycle = std::vector<edge_t> {}; // ???
        while (true)
        {
            auto u = this->_pred[v];
            cycle.push_back(this->_edge[v]); // ???
            v = u;
            if (v == handle)
            {
                break;
            }
        }
        return cycle;
    }

    /*!
     * @brief check if it is really a negative cycle
     *
     * @tparam Container
     * @tparam WeightFn
     * @param handle
     * @param dist
     * @param get_weight
     * @return true
     * @return false
     */
    template <typename Container, typename WeightFn>
    auto __is_negative(const node_t& handle, Container& dist,
        const WeightFn& get_weight) -> bool
    {
        auto v = handle;
        while (true)
        {
            auto u = this->_pred[v];
            auto e = this->_edge[v];
            auto wt = get_weight(this->_G, e); // ???
            if (dist[v] > dist[u] + wt)
            {
                return true;
            }
            v = u;
            if (v == handle)
            {
                break;
            }
        }
        return false;
    }
};
