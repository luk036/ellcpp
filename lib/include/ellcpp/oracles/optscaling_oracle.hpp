// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_OPTSCALING_ORACLE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_OPTSCALING_ORACLE_HPP 1

#include "network_oracle.hpp"
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

using Arr = xt::xarray<double>;

template <typename Graph, typename Edge>
auto constr(Graph& G, Edge& e, const Arr& x) {
    auto u = G.source(e);
    auto v = G.target(e);
    if (u <= v) { // ???
        return x(0) - G[u][v]['cost'];
    }
    return G[u][v]['cost'] - x(1);
}

template <typename Graph, typename Edge>
auto pconstr(Graph& G, Edge& e, const Arr& x) {
    auto u = G.source(e);
    auto v = G.target(e);
    if (u <= v) {
        return Arr{1., 0.};
    }
    return Arr{0., -1.};
}


template <typename Graph, typename Fn_Eval, typename Grad_Fn>
class optscaling_oracle {
  private:
    Graph  &_G;
    network_oracle<Graph,Fn_Eval,Grad_Fn> _network;

  public:
    explicit optscaling_oracle(self, G) :
        _G{G},
        _network{network_oracle<Graph,Fn_Eval,Grad_Fn>(G, constr, pconstr)}
        {}

    auto operator()(const Arr &x, double t) {
        auto [g, f, feasible] = _network(x);
        if (!feasible) {
            return std::tuple{g, f, t};
        }
        auto s = x(0) - x(1);
        auto fj = s - t;
        if (fj < 0.) {
            t = s;
            fj = 0.;
        }
        return std::tuple{Arr{1., -1.}, fj, t};
    }
};

#endif