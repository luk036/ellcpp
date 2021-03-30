// -*- coding: utf-8 -*-
#pragma once

#include <ellcpp/ell.hpp>

/*!
 * @brief Ellipsoid Search Space
 *
 *    ell_stable = {x | (x - xc)' M^-1 (x - xc) \le \kappa}
 *               = {x | (x - xc)' L D^-1 L' (x - xc) \le \kappa}
 *
 * Store $M$ in the form of Lg \ D^-1 \ L' in an n x n array `Q`,
 * and hence keep $M$ symmetric positive definite.
 * More stable but slightly more computation.
 */
class ell_stable : public ell
{
  public:
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

    /*!
     * @brief Construct a new ell_stable object
     *
     * @param[in] val
     * @param[in] x
     */
    ell_stable(const Arr& val, Arr x) noexcept
        : ell {val, std::move(x)}
    {
    }

    /*!
     * @brief Construct a new ell_stable object
     *
     * @param[in] alpha
     * @param[in] x
     */
    ell_stable(const double& alpha, Arr x) noexcept
        : ell {alpha, std::move(x)}
    {
    }

    /**
     * @brief Construct a new ell_stable object
     *
     * @param[in] E (move)
     */
    ell_stable(ell_stable&& E) = default;

  public:
    /**
     * @brief Construct a new ell_stable object
     *
     * @param E
     */
    explicit ell_stable(const ell_stable& E) = default;

    /**
     * @brief explicitly copy
     *
     * @return ell_stable
     */
    [[nodiscard]] auto copy() const -> ell_stable
    {
        return ell_stable(*this);
    }

    /*!
     * @brief Update ellipsoid core function using the cut(s)
     *
     * Overwrite the base class.
     * Store Q^-1 in the form of LDLT decomposition,
     * and hence guarantee Q is symmetric positive definite.
     *
     * @tparam T
     * @param[in] cut cutting-plane
     * @return std::tuple<int, double>
     */
    template <typename T>
    auto update(const std::tuple<Arr, T>& cut) -> std::tuple<CUTStatus, double>;
}; // } ell_stable
