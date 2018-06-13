// -*- coding: utf-8 -*-
/**
Negative cycle detection for (auto weighed graphs.
**/
#include <xnetwork.hpp> // as xn
#include <vector>
#include <unordered_map>

auto default_get_weight(const Graph& G, const Edge& e) {
    auto [u, v] = e;
    return G[u][v].get("weight", 1);
}


class negCycleFinder {
private:
    Graph G;
    std::unordered_map<Node*, int> dist; // int???
    std::unordered_map<Node*, Node*> pred;


public:
    explicit negCycleFinder( const Graph& G, auto get_weight=default_get_weight) {
        this->G = G;
        for (auto v : this->G) {
            this->dist[v] = 0;
            this->pred[v] = nullptr;
        }
        // this->pred = {v: None for v : G};
        this->get_weight = get_weight;
    }

    auto find_cycle( ) {
        /** Find a cycle on policy graph

        Arguments) {
            G {xnetwork graph} 
            pred {dictionary} -- policy graph

        Returns) {
            handle -- a start node of the cycle
        **/

        std::unordered_map<Node*, Node*> visited{};

        for (auto v : this->G) {
            if (visited.find(v) != visited.end()) { // C++20
                continue;
            }
            auto u = v;
            while (true) {
                visited[u] = v;
                u = this->pred[u];
                if (u == nullptr) {
                    break;
                }
                if (visited.find(u) != visited.end()) {
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

    auto relax( ) {
        /** Perform a updating of dist and pred

        Arguments) {
            G {xnetwork graph} -- [description];
            dist {dictionary} -- [description];
            pred {dictionary} -- [description];

        Keyword Arguments) {
            weight {str} -- [description];

        Returns) {
            [type] -- [description];
        **/

        bool changed = false;
        for (auto& e : this->G.edges) {
            auto wt = this->get_weight(this->G, e);
            auto [u, v] = e;
            auto d = this->dist[u] + wt;
            if (this->dist[v] > d) {
                this->dist[v] = d
                this->pred[v] = u;
                changed = true;
            }
        }
        return changed;
    }

    auto find_neg_cycle( ) {
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

    auto neg_cycle_relax( ) {
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

    auto cycle_list( const Node* handle) {
        auto v = handle;
        std::vector<std::tuple<Node*, Node*>> cycle{}; // ???
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

    auto is_negative( const Node* handle) {
        auto v = handle;
        while (true) {
            auto u = this->pred[v];
            auto wt = this->get_weight(this->G, (u, v)); // ???
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
