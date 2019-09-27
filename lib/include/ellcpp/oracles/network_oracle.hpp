// -*- coding: utf-8 -*-
#pragma once

#include "neg_cycle.hpp" // import negCycleFinder
#include <xtensor/xarray.hpp>

/*!
 * @brief
 *
 * @tparam Graph
 * @tparam Fn_Eval
 * @tparam Grad_Fn
 */
template <typename Graph, typename Fn_Eval, typename Grad_Fn>
class network_oracle
{
  private:
    Graph& _G;
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
    explicit network_oracle(Graph& G, Fn_Eval& f, Grad_Fn& p)
        : _G {G}
        , _f {f}
        , _p {p} // partial derivative of f w.r.t x
    {
    }

    /*!
     * @brief
     *
     * @param x
     * @return auto
     */
    auto operator()(const Arr& x) const -> std::tuple<Arr, double>
    {
        auto get_weight = [this, &x](Graph& G, const edge_t& e) -> double {
            return this->_f(G, e, x);
        };

        auto S = negCycleFinder(_G, get_weight);
        auto C = S.find_neg_cycle();
        if (C.empty())
        {
            return {Arr{0.}, -1.};
        }

        auto g = Arr {xt::zeros<double>({x.size()})};
        auto f = 0.;
        for (const auto& e : C)
        {
            f -= _f(_G, e, x);
            g -= _p(_G, e, x);
        }
        return {std::move(g), f};
    }
};

// Template guided deduction
// template <typename Graph, typename Fn_Eval, typename Grad_Fn>
// network_oracle(Graph &G, Fn_Eval &f, Grad_Fn &p)
//     ->network_oracle<Graph, Fn_Eval, Grad_Fn>;
