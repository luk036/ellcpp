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
class network3_oracle
{
  private:
    double _t {0.};

    Graph& _G;
    Fn_Eval _f;
    Grad_Fn _p;

    using edge_t = typename Graph::edge_t;
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

  public:
    /*!
     * @brief Construct a new network oracle object
     *
     * @param G
     * @param f
     * @param p
     */
    explicit network3_oracle(Graph& G, Fn_Eval& f, Grad_Fn& p)
        : _G {G}
        , _f {f}
        , _p {p} // partial derivative of f w.r.t x
    {
    }

    auto update(double t) -> void
    {
        this->_t = t;
    }

    /*!
     * @brief
     *
     * @param x
     * @return auto
     */
    auto operator()(const Arr& x) const
    {
        auto get_weight = [this, &x](Graph& G, const edge_t& e) -> double {
            return this->_f(G, e, x, this->_t);
        };

        auto S = negCycleFinder(_G, get_weight);
        auto C = S.find_neg_cycle();
        if (C.empty())
        {
            return std::tuple {Arr {0.}, -1.};
        }

        auto g = Arr {xt::zeros<double>({x.size()})};
        auto f = 0.;
        for (const auto& e : C)
        {
            f -= this->_f(this->_G, e, x, this->_t);
            g -= this->_p(this->_G, e, x, this->_t);
        }
        return std::tuple {std::move(g), f};
    }
};

// Template guided deduction
// template <typename Graph, typename Fn_Eval, typename Grad_Fn>
// network3_oracle(Graph &G, Fn_Eval &f, Grad_Fn &p)
//     ->network3_oracle<Graph, Fn_Eval, Grad_Fn>;
