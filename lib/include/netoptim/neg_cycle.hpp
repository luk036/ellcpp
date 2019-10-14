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
 * @tparam WeightFn
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
  private:
    const Graph& _G; // const???
    // WeightFn _get_weight; // for nonlinear and lazy evaluation

    using node_t = typename Graph::node_t;
    using edge_t = typename Graph::edge_t;
    // using wt_t = decltype(_get_weight(_G, std::declval<edge_t&>()));

  public:
    py::dict<node_t, node_t> _pred {};
    py::dict<node_t, edge_t> _edge {};
    // py::dict<node_t, wt_t> _dist;

  public:
    /*!
     * @brief Construct a new neg Cycle Finder object
     *
     * @param G
     * @param get_weight
     */
    explicit negCycleFinder(const Graph& G)
        : _G {G}
    {
    }

  private:
    /*!
     * @brief Find a cycle on policy graph
     *
     * @return handle -- a start node of the cycle
     */
    auto find_cycle() -> node_t
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
                        // if (this->is_negative(u)) {
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
     * @brief Perform a updating of dist and pred
     *
     * @return bool
     */
    template <typename Container, typename WeightFn>
    auto relax(Container& dist, const WeightFn& get_weight) -> bool
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

  public:
    /*! Perform a updating of dist and pred
     *    Arguments:
     *        G {[type]} -- [description];
     *        dist {dictionary} -- [description];
     *        pred {dictionary} -- [description];
     *    Keyword Arguments:
     *        weight {str} -- [description] (default: {"weight"});
     *    Returns:
     *        [type] -- [description];
     */
    template <typename Container, typename WeightFn>
    auto find_neg_cycle(Container& dist, const WeightFn& get_weight)
        -> std::vector<edge_t>
    {
        // for (node_t v : this->_G)
        //     this->_dist[v] = wt_t(0);
        this->_pred.clear();
        while (true)
        {
            auto changed = this->relax(dist, get_weight);
            if (!changed)
            {
                break;
            }
            auto v = this->find_cycle();
            if (v != this->_G.null_vertex())
            {
                assert(this->is_negative(v, dist, get_weight));
                return this->cycle_list(v);
            }
        }
        return std::vector<edge_t> {}; // ???
    }

  private:
    /*!
     * @brief
     *
     * @param handle
     * @return std::vector<edge_t>
     */
    auto cycle_list(node_t handle) -> std::vector<edge_t>
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

    template <typename Container, typename WeightFn>
    auto is_negative(const node_t& handle, Container& dist,
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

// Template guided deduction
// template <typename Graph, typename WeightFn>
// negCycleFinder(Graph &G, WeightFn &get_weight)->negCycleFinder<Graph,
// WeightFn>;
