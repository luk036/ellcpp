// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_OPTSCALING_ORACLE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_OPTSCALING_ORACLE_HPP 1

#include "network_oracle.hpp"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

template <typename Graph, typename Fn, typename T> //
class optscaling_oracle {
  private:
    Graph &_G;
    Fn _get_cost;

    using Arr = xt::xarray<double>;
    using edge_t = decltype(*(std::begin(_G.edges())));

  public:
    explicit optscaling_oracle(Graph &G, Fn get_cost, T && /* dummy */)
        : _G{G}, _get_cost{get_cost} {}

    auto operator()(const Arr &x, double t) {
        auto constr = [this](Graph &G, const edge_t &e, const Arr &x) {
            const auto &[u, v] = G.end_points(e);
            if (u <= v) { // ???
                return x(0) - this->_get_cost(G, e);
            }
            return this->_get_cost(G, e) - x(1);
        };

        auto pconstr = [](Graph &G, const edge_t &e, const Arr &) {
            // auto u = G.source(e);
            // auto v = G.target(e);
            const auto &[u, v] = G.end_points(e);
            if (u <= v) {
                return Arr{1., 0.};
            }
            return Arr{0., -1.};
        };

        auto P = network_oracle(_G, constr, pconstr);
        auto [g, f, feasible] = P(x);
        if (!feasible) {
            return std::tuple{std::move(g), f, t};
        }
        auto s = x(0) - x(1);
        auto fj = s - t;
        if (fj < 0) {
            t = s;
            fj = 0.;
        }
        return std::tuple{Arr{1., -1.}, fj, t};
    }
};

#endif