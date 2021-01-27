// -*- coding: utf-8 -*-
#pragma once

#include <ellcpp/utility.hpp>
#include <netoptim/neg_cycle.hpp> // import negCycleFinder
#include <optional>

/*!
 * @brief Oracle for Parametric Network Problem.
 *
 *    This oracle solves the following feasibility problem:
 *
 *        find    x, u
 *        s.t.    u[j] - u[i] \le h(e, x)
 *                \forall e(i, j) \in E
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
    Container& _u; // reference???
    negCycleFinder<Graph> _S;
    Fn _h;

  public:
    /*!
     * @brief Construct a new network oracle object
     *
     * @param[in] G a directed graph (V, E)
     * @param[in,out] u list or dictionary
     * @param[in] h function evaluation and gradient
     */
    network_oracle(const Graph& G, Container& u, Fn h)
        : _G {G}
        , _u {u}
        , _S(G)
        , _h {std::move(h)}
    {
    }

    /**
     * @brief Construct a new network oracle object
     *
     */
    explicit network_oracle(const network_oracle&) = default;

    // network_oracle& operator=(const network_oracle&) = delete;
    // network_oracle(network_oracle&&) = default;

    /*!
     * @brief
     *
     * @param[in] t the best-so-far optimal value
     */
    template <typename opt_type>
    auto update(const opt_type& t) -> void
    {
        this->_h.update(t);
    }

    /*!
     * @brief Make object callable for cutting_plane_feas()
     *
     * @tparam T
     * @param[in] x
     * @return std::optional<std::tuple<T, double>>
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
        for (auto&& e : C)
        {
            f -= this->_h.eval(e, x);
            g -= this->_h.grad(e, x);
        }
        return {{std::move(g), f}};
    }
};
