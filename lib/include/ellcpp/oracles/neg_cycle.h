// -*- coding: utf-8 -*-
/**
Negative cycle detection for (auto weighed graphs.
**/
#include <vector>
// #include <xnetwork.hpp> // as xn
// #include <unordered_map>
#include <py2cpp/py2cpp.hpp>

template <typename Graph>
auto default_get_weight(const Graph &G, const std::tuple<Node *, Node *> &e) {
    auto [u, v] = e;
    return G[u][v].get("weight", 1);
}


template <typename G, typename Fn>
class negCycleFinder {
  private:
    Graph _G;
    Fn _get_weight;
    py::dict<Node *, int> dist; // int???
    py::dict<Node *, Node *> pred;

  public:
    explicit negCycleFinder(const Graph &G,
                            Fn get_weight = default_get_weight)
      : _G{G},
        _get_weight{get_weight} {

        for (auto v : Vertices(_G)) {
            this->dist[v] = 0;
            this->pred[v] = nullptr;
        }
        // this->pred = {v: None for v : G};
    }

    /**
     * @brief Find a cycle on policy graph
     * 
     * @return handle -- a start node of the cycle 
     */
    auto find_cycle() {
        py::dict<Node *, Node *> visited{};

        for (auto v : Vertices(G)) {
            if (visited.contains(v)) {
                continue;
            }
            auto u = v;
            while (true) {
                visited[u] = v;
                u = this->pred[u];
                if (u == nullptr) {
                    break;
                }
                if (visited.contains(u)) {
                    if (visited[u] == v) {
                        if (this->is_negative(u)) {
                            // should be "yield u";
                            return u;
                        }
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
    auto relax() {
        bool changed = false;
        for (auto &e : this->G.edges) {
            auto wt = this->get_weight(this->G, e);
            auto [u, v] = e;
            auto d = this->dist[u] + wt;
            if (this->dist[v] > d) {
                this->dist[v] = d;
                this->pred[v] = u;
                changed = true;
            }
        }
        return changed;
    }

    auto find_neg_cycle() {
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
        for (auto v : this->G) {
            this->dist[v] = 0;
            this->pred[v] = nullptr;
        }
        return this->neg_cycle_relax();
    }

    auto neg_cycle_relax() {
        for (auto v : this->G) {
            this->pred[v] = nullptr;
        }

        while (true) {
            auto changed = this->relax();
            if (changed) {
                // if (v != nullptr) {
                auto v = this->find_cycle();
                if (v != nullptr) {
                    return this->cycle_list(v);
                }
            } else {
                break;
            }
        }
    }

    auto cycle_list(const Node *handle) {
        auto v = handle;
        std::vector<std::tuple<Node *, Node *>> cycle{}; // ???
        while (true) {
            auto u = this->pred[v];
            cycle.emplace_back(std::tuple{u, v});
            v = u;
            if (v == handle) {
                break;
            }
        }
        return cycle;
    }

    auto is_negative(const Node *handle) {
        auto v = handle;
        while (true) {
            auto u = this->pred[v];
            auto wt = this->get_weight(this->G, {u, v}); // ???
            if (this->dist[v] > this->dist[u] + wt) {
                return true;
            }
            v = u;
            if (v == handle) {
                break;
            }
        }
        return false;
    }
