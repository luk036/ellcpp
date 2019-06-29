// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_NEG_CYCLE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_NEG_CYCLE_HPP 1

/*!
Negative cycle detection for (auto weighed graphs.
**/
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
template<typename Graph, typename WeightFn>
class negCycleFinder
{
private:
    Graph&   _G;
    WeightFn _get_weight; // for nonlinear and lazy evaluation

    using Node   = decltype(Graph::null_vertex());
    using edge_t = decltype(*(std::begin(_G.edges())));
    using wt_t   = decltype(_get_weight(_G, std::declval<edge_t&>()));

public:
    py::dict<Node, Node>   _pred;
    py::dict<Node, edge_t> _edge;
    py::dict<Node, wt_t>   _dist;

public:
    /*!
     * @brief Construct a new neg Cycle Finder object
     *
     * @param G
     * @param get_weight
     */
    explicit negCycleFinder(Graph& G, WeightFn& get_weight) : _G{G}, _get_weight{get_weight}
    {
        for (Node v : _G)
        {
            _dist[v] = wt_t(0);
        }
        _pred.clear();
        // _pred = {v: None for v : _G};
    }

    /*!
     * @brief Find a cycle on policy graph
     *
     * @return handle -- a start node of the cycle
     */
    Node find_cycle()
    {
        py::dict<Node, Node> visited{};

        for (Node v : _G)
        {
            if (visited.contains(v)) { continue; }
            auto u = v;
            while (true)
            {
                visited[u] = v;
                if (!_pred.contains(u)) { break; }
                u = _pred[u];
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

        return _G.null_vertex();
    }

    /*!
     * @brief Perform a updating of dist and pred
     *
     * @return auto
     */
    auto relax() -> bool
    {
        auto changed = false;
        for (auto const& e : _G.edges())
        {
            auto wt       = _get_weight(_G, e);
            auto&& [u, v] = _G.end_points(e);
            auto d        = _dist[u] + wt;
            if (_dist[v] > d)
            {
                _dist[v] = d;
                _pred[v] = u;
                _edge[v] = e;
                changed  = true;
            }
        }
        return changed;
    }

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
    auto find_neg_cycle() -> std::vector<edge_t>
    {
        for (Node v : _G)
        {
            _dist[v] = wt_t(0);
        }
        _pred.clear();
        return std::move(this->neg_cycle_relax());
    }

    /*!
     * @brief
     *
     * @return std::vector<edge_t>
     */
    auto neg_cycle_relax() -> std::vector<edge_t>
    {
        _pred.clear();

        while (true)
        {
            auto changed = this->relax();
            if (!changed) { break; }
            // if (v != _G.null_vertex()) {
            Node v = this->find_cycle();
            if (v != _G.null_vertex()) { return std::move(this->cycle_list(v)); }
        }
        return std::move(std::vector<edge_t>{}); // ???
    }

    /*!
     * @brief
     *
     * @param handle
     * @return std::vector<edge_t>
     */
    auto cycle_list(Node handle) -> std::vector<edge_t>
    {
        auto v     = handle;
        auto cycle = std::vector<edge_t>{}; // ???
        while (true)
        {
            auto u = _pred[v];
            cycle.push_back(_edge[v]);
            v = u;
            if (v == handle) { break; }
        }
        return std::move(cycle);
    }

    // auto is_negative(const Node &handle) {
    //     Node v = handle;
    //     while (true) {
    //         auto u = _pred[v];
    //         auto wt = _get_weight(_G, {u, v}); // ???
    //         if (_dist[v] > _dist[u] + wt) {
    //             return true;
    //         }
    //         v = u;
    //         if (v == handle) {
    //             break;
    //         }
    //     }
    //     return false;
    // }
};

// Template guided deduction
// template <typename Graph, typename WeightFn>
// negCycleFinder(Graph &G, WeightFn &get_weight)->negCycleFinder<Graph,
// WeightFn>;

#endif