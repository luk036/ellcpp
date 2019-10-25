// -*- coding: utf-8 -*-
#pragma once

#include "network_oracle.hpp"
#include <cassert>
#include <xtensor/xarray.hpp>

/**
 * @brief 
 * 
 * @tparam Graph 
 * @tparam Container 
 * @tparam Fn 
 */
template <typename Graph, typename Container, typename Fn> //
class optscaling_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using edge_t = typename Graph::edge_t;
    using Cut = std::tuple<Arr, double>;

    // using constr_fn = std::function<double(const Graph&, const edge_t, const
    // Arr&)>; using pconstr_fn = std::function<Arr(const Graph&, const edge_t,
    // const Arr&)>; using NWO = network_oracle<Graph, Container, constr_fn,
    // pconstr_fn>;

  public:
    struct Ratio
    {
        Fn _get_cost;

        explicit Ratio(Fn get_cost)
            : _get_cost {get_cost}
        {
        }

        auto eval(const Graph& G, const edge_t& e, const Arr& x) const
            -> double
        {
            auto [u, v] = G.end_points(e);
            auto cost = this->_get_cost(G, e);
            assert(u != v);
            return (u < v) ? x(0) - cost : cost - x(1);
        }

        auto grad(const Graph& G, const edge_t& e, const Arr& x) const
            -> Arr
        {
            auto [u, v] = G.end_points(e);
            assert(u != v);
            return (u < v) ? Arr {1., 0.} : Arr {0., -1.};
        }
    };

  private:
    network_oracle<Graph, Container, Ratio> _network;

  public:
    /*!
     * @brief Construct a new optscaling oracle object
     *
     * @param G
     * @param get_cost
     */
    optscaling_oracle(const Graph& G, Container& dist, Fn get_cost)
        : _network(G, dist, Ratio {get_cost})
    {
    }

    /*!
     * @brief
     *
     * @param x
     * @param t
     * @return auto
     */
    auto operator()(const Arr& x, double t) -> std::tuple<Cut, double>
    {
        if (auto cut = this->_network(x))
        {
            return {*cut, t};
        }
        auto s = x(0) - x(1);
        auto fj = s - t;
        if (fj < 0)
        {
            t = s;
            fj = 0.;
        }
        return {{Arr {1., -1.}, fj}, t};
    }
};
