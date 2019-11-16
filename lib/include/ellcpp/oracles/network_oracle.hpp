// -*- coding: utf-8 -*-
#pragma once

#include <ellcpp/utility.hpp>
#include <netoptim/neg_cycle.hpp> // import negCycleFinder
#include <optional>

/*!
 * @brief Oracle for Parametric Network Problem:
 *
 *        find    x, u
 *        s.t.    u[j] - u[i] ≤ h(e, x)
 *                ∀ e(i, j) ∈ E
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
    Container& _u;
    Fn _h;
    negCycleFinder<Graph> _S;

  public:
    /*!
     * @brief Construct a new network oracle object
     *
     * @param G a directed graph (V, E)
     * @param u list or dictionary
     * @param h function evaluation and gradient
     */
    network_oracle(const Graph& G, Container& u, const Fn& h)
        : _G {G}
        , _u {u}
        , _h {h}
        , _S(G)
    {
    }

    /*!
     * @brief
     *
     * @param t the best-so-far optimal value
     * @return auto
     */
    auto update(const double& t)
    {
        this->_h.update(t);
    }

    /*!
     * @brief Make object callable for cutting_plane_feas()
     *
     * @param x
     * @return auto
     */
    template <typename T>
    auto operator()(const T& x) -> std::optional<std::tuple<T, double>>
    {
        auto get_weight = [this, &x](const edge_t& e) -> double {
            return this->_h.eval(e, x);
        };

        auto C = this->_S.find_neg_cycle(this->_u, get_weight);
        if (C.empty())
        {
            return {};
        }

        auto g = zeros(x);
        auto f = 0.;
        for (const auto& e : C)
        {
            f -= this->_h.eval(e, x);
            g -= this->_h.grad(e, x);
        }
        return {{std::move(g), f}};
    }
};
