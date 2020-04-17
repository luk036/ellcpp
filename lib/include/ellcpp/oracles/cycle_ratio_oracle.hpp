// -*- coding: utf-8 -*-
#pragma once

#include "network_oracle.hpp"
#include <cassert>

/*!
 * @brief Oracle for minimum ratio cycle problem.
 *
 *    This example solves the following convex problem:
 *
 *        max     t
 *        s.t.    u[j] - u[i] \le mij - sij * x,
 *                t \le x
 *
 *    where sij is not necessarily positive.
 *
 * @tparam Graph
 * @tparam Container
 * @tparam Fn1
 * @tparam Fn2
 */
template <typename Graph, typename Container, typename Fn1, typename Fn2> //
class cycle_ratio_oracle
{
    using edge_t = typename Graph::edge_t;
    using Cut = std::tuple<double, double>;

    /**
     * @brief Ratio
     *
     */
    class Ratio
    {
      private:
        const Graph& _G;
        Fn1 _get_cost;
        Fn2 _get_time;

      public:
        /*!
         * @brief Construct a new Ratio object
         *
         * @param[in] G
         * @param[in] get_cost
         * @param[in] get_time
         */
        Ratio(const Graph& G, Fn1 get_cost, Fn2 get_time)
            : _G {G}
            , _get_cost {std::move(get_cost)}
            , _get_time {std::move(get_time)}
        {
        }

        /*!
         * @brief Evaluate function
         *
         * @param[in] e
         * @param[in] x (\pi, \phi) in log scale
         * @return double
         */
        auto eval(const edge_t& e, const double& x) const -> double
        {
            const auto cost = this->_get_cost(e);
            const auto time = this->_get_time(e);
            return cost - time * x;
        }

        /*!
         * @brief Gradient function
         *
         * @param[in] e
         * @param[in] x (\pi, \phi) in log scale
         * @return double
         */
        auto grad(const edge_t& e, const double& x) const -> double
        {
            return -this->_get_time(e);
        }
    };

    network_oracle<Graph, Container, Ratio> _network;

  public:
    /*!
     * @brief Construct a new cycle_ratio oracle object
     *
     * @param[in] G
     * @param[in,out] u
     * @param[in] get_cost
     */
    cycle_ratio_oracle(const Graph& G, Container& u, Fn1 get_cost, Fn2 get_time)
        : _network(G, u, Ratio {G, get_cost, get_time})
    {
    }

    cycle_ratio_oracle(const cycle_ratio_oracle&) = delete;
    cycle_ratio_oracle& operator=(const cycle_ratio_oracle&) = delete;
    cycle_ratio_oracle(cycle_ratio_oracle&&) = default;

    /*!
     * @brief Make object callable for cutting_plane_dc()
     *
     * @param[in] x
     * @param[in] t the best-so-far optimal value
     * @return std::tuple<Cut, double>
     *
     * @see cutting_plane_dc
     */
    auto operator()(const double& x, double t) -> std::tuple<Cut, double>
    {
        auto fj = t - x;
        if (fj >= 0)
        {
            return {{-1., fj}, t};
        }
        if (auto cut = this->_network(x))
        {
            return {*cut, t};
        }
        return {{-1., 0.}, x};
    }
};
