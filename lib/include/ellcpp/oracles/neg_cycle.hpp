// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_NEG_CYCLE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_NEG_CYCLE_HPP 1

/**
Negative cycle detection for (auto weighed graphs.
**/
#include <vector>
#include <py2cpp/py2cpp.hpp>

/**
 * @brief negative cycle
 *
 * @tparam Graph
 * @tparam WeightFn
 */
template <typename Graph, typename WeightFn, typename wt_t>
class negCycleFinder {
    using Node = decltype(Graph::null_vertex());
    using edge_t = decltype(*(std::declval<Graph>().edges().begin()));
    // using wt_t = decltype( boost::get(boost::edge_weight,
    // std::declval<Graph>() )[std::declval<edge_t>()] );

  private:
    Graph &_G;
    WeightFn _get_weight;

  public:
    py::dict<Node, wt_t>    _dist;
    py::dict<Node, Node>    _pred;
    py::dict<Node, edge_t>  _edge;

  public:
    explicit negCycleFinder(Graph &G, WeightFn &get_weight,
                            wt_t /* better idea??? */) :
        _G{G}, 
        _get_weight{get_weight} {

        for (Node v : _G) {
            _dist[v] = wt_t(0);
            _pred[v] = _G.null_vertex();
        }
        // _pred = {v: None for v : _G};
    }

    /**
     * @brief Find a cycle on policy graph
     *
     * @return handle -- a start node of the cycle
     */
    Node find_cycle() {
        py::dict<Node, Node> visited{};

        for (Node v : _G) {
            if (visited.contains(v)) {
                continue;
            }
            auto u = v;
            while (true) {
                visited[u] = v;
                u = _pred[u];
                if (u == _G.null_vertex()) {
                    break;
                }
                if (visited.contains(u)) {
                    if (visited[u] == v) {
                        // if (this->is_negative(u)) {
                        // should be "yield u";
                        return u;
                        // }
                    }
                    break;
                }
            }
        }
    }

    /**
     * @brief Perform a updating of dist and pred
     *
     * @return auto
     */
    bool relax() {
        bool changed = false;
        for (auto e : _G.edges()) {
            auto wt = _get_weight(_G, e);
            // auto [u, v] = e;
            auto u = _G.source(e);
            auto v = _G.target(e);
            auto d = _dist[u] + wt;
            if (_dist[v] > d) {
                _dist[v] = d;
                _pred[v] = u;
                _edge[v] = e;
                changed = true;
            }
        }
        return changed;
    }

    /** Perform a updating of dist and pred
     *    Arguments:
     *        G {[type]} -- [description];
     *        dist {dictionary} -- [description];
     *        pred {dictionary} -- [description];
     *    Keyword Arguments:
     *        weight {str} -- [description] (default: {"weight"});
     *    Returns:
     *        [type] -- [description];
     */
    auto find_neg_cycle() {
        for (Node v : _G) {
            _dist[v] = wt_t(0);
            _pred[v] = _G.null_vertex();
        }
        return this->neg_cycle_relax();
    }

    auto neg_cycle_relax() {
        for (Node v : _G) {
            _pred[v] = _G.null_vertex();
        }

        while (true) {
            auto changed = this->relax();
            if (!changed) {
                break;
            }
            // if (v != _G.null_vertex()) {
            Node v = this->find_cycle();
            if (v != _G.null_vertex()) {
                return this->cycle_list(v);
            }
        }
        return std::vector<edge_t>{}; // ???
    }

    auto cycle_list(Node handle) {
        Node v = handle;
        std::vector<edge_t> cycle{}; // ???
        while (true) {
            auto u = _pred[v];
            cycle.push_back(_edge[v]);
            v = u;
            if (v == handle) {
                break;
            }
        }
        return cycle;
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

#endif