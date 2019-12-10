// -*- coding: utf-8 -*-
#pragma once

#include "network_oracle.hpp"
#include <cassert>
#include <xtensor/xarray.hpp>

/*!
 * @brief Oracle for Optimal Matrix Scaling.
 *
 *    This example is taken from [Orlin and Rothblum, 1985]:
 *
 *        min     π/ψ
 *        s.t.    ψ ≤ u[i] * |aij| * u[j]^−1 ≤ π,
 *                ∀ aij != 0,
 *                π, ψ, u, positive
 *
 * @tparam Graph
 * @tparam Container
 * @tparam Fn
 */
template <typename Graph, typename Container, typename Fn> //
class optscaling_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using edge_t = typename Graph::edge_t;
    using Cut = std::tuple<Arr, double>;

  public:
    /**
     * @brief Ratio
     *
     */
    struct Ratio
    {
        const Graph& _G;
        Fn _get_cost;

        /*!
         * @brief Construct a new Ratio object
         *
         * @param G
         * @param get_cost
         */
        Ratio(const Graph& G, Fn get_cost)
            : _G {G}
            , _get_cost {std::move(get_cost)}
        {
        }

        Ratio& operator=(const Ratio& ) = delete;

        /*!
         * @brief Evaluate function
         *
         * @param e
         * @param x (π, ψ) in log scale
         * @return double
         */
        auto eval(const edge_t& e, const Arr& x) const -> double
        {
            const auto [u, v] = this->_G.end_points(e);
            const auto cost = this->_get_cost(e);
            assert(u != v);
            return (u < v) ? x(0) - cost : cost - x(1);
        }

        /*!
         * @brief Gradient function
         *
         * @param e
         * @param x (π, ψ) in log scale
         * @return Arr
         */
        auto grad(const edge_t& e, const Arr& x) const -> Arr
        {
            const auto [u, v] = this->_G.end_points(e);
            assert(u != v);
            return (u < v) ? Arr {1., 0.} : Arr {0., -1.};
        }
    };

  private:
    network_oracle<Graph, Container, Ratio> _network;

  public:
    /*!
     * @brief Construct a new optscaling oracle object
     *
     * @param G
     * @param u
     * @param get_cost
     */
    optscaling_oracle(const Graph& G, Container& u, Fn get_cost)
        : _network(G, u, Ratio {G, get_cost})
    {
    }

    optscaling_oracle(const optscaling_oracle& ) = delete;
    optscaling_oracle& operator=(const optscaling_oracle& ) = delete;

    /*!
     * @brief Make object callable for cutting_plane_dc()
     *
     * @param x (π, ψ) in log scale
     * @param t the best-so-far optimal value
     * @return std::tuple<Cut, double>
     *
     * @see cutting_plane_dc
     */
    auto operator()(const Arr& x, double t) -> std::tuple<Cut, double>
    {
        if (auto cut = this->_network(x))
        {
            return {*cut, t};
        }
        auto s = x(0) - x(1);
        auto fj = s - t;
        if (fj < 0)
        {
            t = s;
            fj = 0.;
        }
        return {{Arr {1., -1.}, fj}, t};
    }
};
