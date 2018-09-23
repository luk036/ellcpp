// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_NETWORK_ORACLE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_NETWORK_ORACLE_HPP 1

#include "neg_cycle.hpp" // import negCycleFinder
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

template <typename Graph, typename Fn_Eval, typename Grad_Fn>
class network_oracle {
  private:
    Graph    &_G;
    Fn_Eval  _f;
    Grad_Fn  _p;

    using edge_t = decltype(*(_G.edges().begin()));

  public:
    explicit network_oracle(Graph &G, Fn_Eval &f, Grad_Fn &p) :
        _G{G},
        _f{f},
        _p{p}  // partial derivative of f w.r.t x
        {}

    auto operator()(const Arr &x) {
        auto get_weight = [](Graph &G, edge_t &e) -> double {
            return _f(G, e, x);
        };

        auto S = negCycleFinder(G, get_weight);
        auto C = S.find_neg_cycle();

        auto n = x.size();
        Arr g = xt::zeros<double>({n});
        double f = 0.0;

        if (C.empty()) {
            return std::tuple{g, 0., true};
        }
        for (auto e : C) {
            f -= _f(_G, e, x);
            g -= _p(_G, e, x);
        }
        return std::tuple{g, f, false};
    }
};

#endif