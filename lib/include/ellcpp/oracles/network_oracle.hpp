// -*- coding: utf-8 -*-
#pragma once

#include <netoptim/neg_cycle.hpp> // import negCycleFinder
// #include <xtensor/xarray.hpp>
#include <ellcpp/utility.hpp>

/*!
 * @brief
 *
 * @tparam Graph
 * @tparam Container
 * @tparam Fn_Eval
 * @tparam Grad_Fn
 */
template <typename Graph, typename Container, typename Fn_Eval,
    typename Grad_Fn>
class network_oracle
{
    // using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using edge_t = typename Graph::edge_t;

  private:
    const Graph& _G;
    Container& _dist;
    Fn_Eval _f;
    Grad_Fn _p;
    negCycleFinder<Graph> _S;

  public:
    /*!
     * @brief Construct a new network oracle object
     *
     * @param G
     * @param f
     * @param p
     */
    network_oracle(
        const Graph& G, Container& dist, const Fn_Eval& f, const Grad_Fn& p)
        : _G {G}
        , _dist {dist}
        , _f {f}
        , _p {p} // partial derivative of f w.r.t x
        , _S(G)
    {
    }

    /*!
     * @brief
     *
     * @param x
     * @return auto
     */
    template <typename T>
    auto operator()(const T& x)
    {
        auto get_weight = [this, &x](
                              const Graph& G, const edge_t& e) -> double {
            return this->_f(G, e, x);
        };

        // auto S = negCycleFinder(this->_G);
        auto C = this->_S.find_neg_cycle(this->_dist, get_weight);
        if (C.empty())
        {
            return std::tuple {false, T {0.}, -1.};
        }

        auto g = zeros(x);
        auto f = 0.;
        for (const auto& e : C)
        {
            f -= this->_f(this->_G, e, x);
            g -= this->_p(this->_G, e, x);
        }
        return std::tuple {true, std::move(g), f};
    }
};

// Template guided deduction
// template <typename Graph, typename Fn_Eval, typename Grad_Fn>
// network_oracle(Graph &G, Fn_Eval &f, Grad_Fn &p)
//     ->network_oracle<Graph, Fn_Eval, Grad_Fn>;
