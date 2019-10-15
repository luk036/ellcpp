// -*- coding: utf-8 -*-
#pragma once

#include "network_oracle.hpp"
#include <xtensor/xarray.hpp>
#include <functional>
#include <memory>

/*!
 * @brief
 *
 * @tparam Graph
 * @tparam Fn
 * @tparam T
 */
template <typename Graph, typename Container, typename Fn> //
class optscaling_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using edge_t = typename Graph::edge_t;
    // using constr_fn = std::function<double(const Graph&, const edge_t, const Arr&)>;
    // using pconstr_fn = std::function<Arr(const Graph&, const edge_t, const Arr&)>;
    // using NWO = network_oracle<Graph, Container, constr_fn, pconstr_fn>;

  public:
    struct constr
    {
        Fn _get_cost;

        explicit constr(Fn get_cost)
            : _get_cost {get_cost}
        {
        }

        auto operator()(const Graph& G, const edge_t& e, const Arr& x) const -> double
        {
            auto [u, v] = G.end_points(e);
            auto cost = this->_get_cost(G, e);
            return (u <= v) ? x(0) - cost : cost - x(1);
        }
    };

    struct pconstr
    {
        auto operator()(const Graph& G, const edge_t& e, const Arr& x) const -> Arr
        {
            auto [u, v] = G.end_points(e);
            return (u <= v) ? Arr {1., 0.} : Arr {0., -1.};
        }
    };

  private:
    // const Graph& _G;
    // Container& _dist;
    // Fn _get_cost;
    network_oracle<Graph, Container, constr, pconstr> _network;

  public:
    /*!
     * @brief Construct a new optscaling oracle object
     *
     * @param G
     * @param get_cost
     */
    optscaling_oracle(const Graph& G, Container& dist, Fn get_cost)
        : _network(G, dist, constr{get_cost}, pconstr{})
    {
    }

    /*!
     * @brief
     *
     * @param x
     * @param t
     * @return auto
     */
    auto operator()(const Arr& x, double t)
    {
        // auto constr = [this](const Graph& G, const edge_t& e, const Arr& x) {
        //     auto [u, v] = G.end_points(e);
        //     auto cost = this->_get_cost(G, e);
        //     return (u <= v) ? x(0) - cost : cost - x(1);
        // };

        // auto pconstr = [](const Graph& G, const edge_t& e, const Arr&) {
        //     auto [u, v] = G.end_points(e);
        //     return (u <= v) ? Arr {1., 0.} : Arr {0., -1.};
        // };

        // auto network = network_oracle(this->_G, this->_dist, constr{this->_get_cost}, pconstr{});
        auto [g, f] = this->_network(x);
        if (g.shape()[0] > 1 || g(0) != 0)
        {
            return std::tuple {std::move(g), f, t};
        }
        auto s = x(0) - x(1);
        auto fj = s - t;
        if (fj < 0)
        {
            t = s;
            fj = 0.;
        }
        return std::tuple {Arr {1., -1.}, fj, t};
    }
};
