// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_NETWORK_ORACLE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ORACLES_NETWORK_ORACLE_HPP 1

#include "neg_cycle.hpp" // import negCycleFinder
// #include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

/*!
 * @brief
 *
 * @tparam Graph
 * @tparam Fn_Eval
 * @tparam Grad_Fn
 */
template <typename Graph, typename Fn_Eval, typename Grad_Fn>
class network_oracle {
  private:
    Graph &_G;
    Fn_Eval _f;
    Grad_Fn _p;

    using edge_t = decltype(*(std::begin(_G.edges())));
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

  public:
    /*!
     * @brief Construct a new network oracle object
     *
     * @param G
     * @param f
     * @param p
     */
    explicit network_oracle(Graph &G, Fn_Eval &f, Grad_Fn &p)
        : _G{G}, _f{f}, _p{p} // partial derivative of f w.r.t x
    {}

    /*!
     * @brief
     *
     * @param x
     * @return auto
     */
    auto operator()(const Arr &x) const {
        auto get_weight = [this, &x](Graph &G, const edge_t &e) -> double {
            return this->_f(G, e, x);
        };

        auto S = negCycleFinder(_G, get_weight);
        auto C = S.find_neg_cycle();
        auto g = Arr{xt::zeros<double>({x.size()})};
        auto f = 0.;

        if (C.empty()) {
            return std::tuple{std::move(g), 0., true};
        }
        for (auto const &e : C) {
            f -= _f(_G, e, x);
            g -= _p(_G, e, x);
        }
        return std::tuple{std::move(g), f, false};
    }
};

// Template guided deduction
// template <typename Graph, typename Fn_Eval, typename Grad_Fn>
// network_oracle(Graph &G, Fn_Eval &f, Grad_Fn &p)
//     ->network_oracle<Graph, Fn_Eval, Grad_Fn>;

#endif
