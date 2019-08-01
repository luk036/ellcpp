// -*- coding: utf-8 -*-
#pragma once

#include "network_oracle.hpp"
#include <xtensor/xarray.hpp>

/*!
 * @brief
 *
 * @tparam Graph
 * @tparam Fn
 * @tparam T
 */
template<typename Graph, typename Fn, typename T> //
class optscaling_oracle
{
private:
    Graph& _G;
    Fn     _get_cost;

    using Arr    = xt::xarray<double, xt::layout_type::row_major>;
    using edge_t = decltype(*(std::begin(_G.edges())));

public:
    /*!
     * @brief Construct a new optscaling oracle object
     *
     * @param G
     * @param get_cost
     */
    explicit optscaling_oracle(Graph& G, Fn get_cost, T&& /* dummy */) : _G{G}, _get_cost{get_cost}
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
        auto constr = [this](Graph& G, const edge_t& e, const Arr& x) {
            auto&& [u, v] = G.end_points(e);
            auto cost     = this->_get_cost(G, e);
            return (u <= v) ? x(0) - cost : cost - x(1);
        };

        auto pconstr = [](Graph& G, const edge_t& e, const Arr&) {
            auto&& [u, v] = G.end_points(e);
            return (u <= v) ? Arr{1., 0.} : Arr{0., -1.};
        };

        auto P                = network_oracle(_G, constr, pconstr);
        auto [g, f, feasible] = P(x);
        if (!feasible) { return std::tuple{std::move(g), f, t}; }
        auto s  = x(0) - x(1);
        auto fj = s - t;
        if (fj < 0)
        {
            t  = s;
            fj = 0.;
        }
        return std::tuple{Arr{1., -1.}, fj, t};
    }
};
