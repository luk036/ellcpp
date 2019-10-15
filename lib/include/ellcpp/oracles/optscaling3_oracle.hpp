// -*- coding: utf-8 -*-
#pragma once

#include "network3_oracle.hpp"
#include <cassert>
// #include <xtensor/xarray.hpp>

/*!
 * @brief
 *
 * @tparam Graph
 * @tparam Fn
 * @tparam T
 */
template <typename Graph, typename Container, typename Fn> //
class optscaling_oracle3
{
    // using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using edge_t = typename Graph::edge_t;
    // using constr_fn = std::function<double(const Graph&, const edge_t, const
    // Arr&)>; using pconstr_fn = std::function<Arr(const Graph&, const edge_t,
    // const Arr&)>; using NWO = network_oracle<Graph, Container, constr_fn,
    // pconstr_fn>;

  public:
    struct constr3
    {
        Fn _get_cost;

        explicit constr3(Fn get_cost)
            : _get_cost {get_cost}
        {
        }

        auto operator()(const Graph& G, const edge_t& e, const double& x,
            double t) const -> double
        {
            auto [u, v] = G.end_points(e);
            auto cost = this->_get_cost(G, e);
            assert(u != v);
            return (u < v) ? x + t - cost : cost - x;
        }
    };

    struct pconstr3
    {
        auto operator()(const Graph& G, const edge_t& e, const double& /* x */,
            double /* t */) const -> double
        {
            auto [u, v] = G.end_points(e);
            assert(u != v);
            return (u < v) ? 1. : -1.;
        }
    };

  private:
    // const Graph& _G;
    // Container& _dist;
    // Fn _get_cost;
    network3_oracle<Graph, Container, constr3, pconstr3> _network3;

  public:
    /*!
     * @brief Construct a new optscaling oracle object
     *
     * @param G
     * @param get_cost
     */
    optscaling_oracle3(const Graph& G, Container& dist, Fn get_cost)
        : _network3(G, dist, constr3 {get_cost}, pconstr3 {})
    {
    }

    auto update(double t) -> void
    {
        this->_network3.update(t);
    }

    /*!
     * @brief
     *
     * @param x
     * @param t
     * @return auto
     */
    auto operator()(const double& x)
    {
        return this->_network3(x);
    }
};
