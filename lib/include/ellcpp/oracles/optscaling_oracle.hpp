// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_OPTSCALING_ORACLE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_OPTSCALING_ORACLE_HPP 1

#include "network_oracle.hpp"
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

template <typename Graph, typename CostFn> class optscaling_oracle {
  private:
    Graph &_G;
    CostFn _get_cost;

    using Arr = xt::xarray<double>;
    using edge_t = decltype(*(_G.edges().begin()));

  public:
    explicit optscaling_oracle(Graph &G, CostFn get_cost)
        : _G{G}, _get_cost{get_cost} {}

    auto operator()(const Arr &x, double t) {
        auto constr = [this](Graph &G, edge_t &e, const Arr &x) {
            auto u = G.source(e);
            auto v = G.target(e);
            if (u <= v) { // ???
                return x(0) - this->_get_cost(G, e, x);
            }
            return this->_get_cost(G, e, x) - x(1);
        };

        auto pconstr = [](Graph &G, edge_t &e, const Arr &) {
            auto u = G.source(e);
            auto v = G.target(e);
            if (u <= v) {
                return Arr{1., 0.};
            }
            return Arr{0., -1.};
        };

        auto P = network_oracle(_G, constr, pconstr);
        auto [g, f, feasible] = P(x);
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