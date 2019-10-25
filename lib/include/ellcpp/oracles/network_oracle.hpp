// -*- coding: utf-8 -*-
#pragma once

#include <ellcpp/utility.hpp>
#include <netoptim/neg_cycle.hpp> // import negCycleFinder
#include <optional>
// #include <xtensor/xarray.hpp>

/*!
 * @brief
 *
 * @tparam Graph
 * @tparam Container
 * @tparam Fn
 */
template <typename Graph, typename Container, typename Fn>
class network_oracle
{
    using edge_t = typename Graph::edge_t;
    
  private:
    const Graph& _G;
    Container& _dist;
    Fn _h;
    negCycleFinder<Graph> _S;

  public:
    /*!
     * @brief Construct a new network oracle object
     *
     * @param G
     * @param f
     * @param p
     */
    network_oracle(const Graph& G, Container& dist, const Fn& h)
        : _G {G}
        , _dist {dist}
        , _h {h}
        , _S(G)
    {
    }

    auto update(const double& t)
    {
        this->_h.update(t);
    }

    /*!
     * @brief
     *
     * @param x
     * @return auto
     */
    template <typename T>
    auto operator()(const T& x) -> std::optional<std::tuple<T, double>>
    {
        auto get_weight = [this, &x](
                              const Graph& G, const edge_t& e) -> double {
            return this->_h.eval(G, e, x);
        };

        // auto S = negCycleFinder(this->_G);
        auto C = this->_S.find_neg_cycle(this->_dist, get_weight);
        if (C.empty())
        {
            return {};
        }

        auto g = zeros(x);
        auto f = 0.;
        for (const auto& e : C)
        {
            f -= this->_h.eval(this->_G, e, x);
            g -= this->_h.grad(this->_G, e, x);
        }
        return {{std::move(g), f}};
    }
};

// Template guided deduction
// template <typename Graph, typename Fn_Eval, typename Grad_Fn>
// network_oracle(Graph &G, Fn_Eval &f, Grad_Fn &p)
//     ->network_oracle<Graph, Fn_Eval, Grad_Fn>;
