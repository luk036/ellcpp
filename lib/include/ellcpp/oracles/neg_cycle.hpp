// -*- coding: utf-8 -*-
/**
Negative cycle detection for (auto weighed graphs.
**/
#include <vector>
// #include <xnetwork.hpp> // as xn
// #include <unordered_map>
#include <py2cpp/py2cpp.hpp>
#include <py2cpp/nx2bgl.hpp>

template <typename Graph, typename WeightFn>
class negCycleFinder {
    using Node = decltype(Graph::null_vertex());
    using Edge = std::pair<Node, Node>;

  private:
    Graph& _G;
    WeightFn _get_weight;
    py::dict<Node, int> _dist; // int???
    py::dict<Node, Node> _pred;
    //EdgeWeightPMap _weightmap;

  public:
    explicit negCycleFinder(Graph &G, WeightFn& get_weight)
        : _G{G}, 
          _get_weight{ get_weight } {

        for (Node v : _G) {
            _dist[v] = 0;
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
        for (auto e : _G.edges() ) {
            int wt = this->_get_weight(_G, e);
            // auto [u, v] = e;
            Node u = _G.source(e);
            Node v = _G.target(e);
            auto d = _dist[u] + wt;
            if (_dist[v] > d) {
                _dist[v] = d;
                _pred[v] = u;
                changed = true;
            }
        }
        return changed;
    }

    std::vector<Edge> find_neg_cycle() {
        /** Perform a updating of dist and pred

        Arguments) {
            G {[type]} -- [description];
            dist {dictionary} -- [description];
            pred {dictionary} -- [description];

        Keyword Arguments) {
            weight {str} -- [description] (default: {"weight"});

        Returns) {
            [type] -- [description];
        **/
        for (Node v : _G) {
            _dist[v] = 0;
            _pred[v] = _G.null_vertex();
        }
        return this->neg_cycle_relax();
    }

    std::vector<Edge> neg_cycle_relax() {
        for (Node v : _G) {
            _pred[v] = _G.null_vertex();
        }

        while (true) {
            auto changed = this->relax();
            if (changed) {
                // if (v != _G.null_vertex()) {
                Node v = this->find_cycle();
                if (v != _G.null_vertex()) {
                    return this->cycle_list(v);
                }
            } else {
                break;
            }
        }
        return std::vector<Edge>{}; // ???

    }

    auto cycle_list(Node handle) {
        Node v = handle;
        std::vector<Edge> cycle{}; // ???
        while (true) {
            auto u = _pred[v];
            cycle.emplace_back(Edge{u, v});
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