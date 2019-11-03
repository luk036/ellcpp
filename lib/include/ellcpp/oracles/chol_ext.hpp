// -*- coding: utf-8 -*-
#pragma once

#include <ellcpp/utility.hpp>
#include <stdexcept>
#include <xtensor/xarray.hpp>

/*!
 * @brief Cholesky factorization for LMI
 *
 *  - LDLT square-root-free version
 *  - Option allow semidefinite
 *  - A matrix $A in R^{m x m}$ is positive definite iff v' A v > 0
 *      for all v in R^n.
 *  - O($p^2 n$) per iteration, independent of $m$
 */
template <bool Allow_semidefinite = false> //
class chol_ext
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Vec = Arr;
    using Mat = Arr;
    using Rng = std::pair<size_t, size_t>;

  public:
    Rng p {0, 0}; /**< the rows where the process starts and stops */
    Vec v;        /**< witness vector */

  private:
    size_t n; /**< dimension */
    Mat T;    /**< temporary storage */

  public:
    /*!
     * @brief Construct a new chol ext object
     *
     * @param N dimension
     */
    explicit chol_ext(size_t N)
        : v {zeros({N})}
        , n {N}
        , T {zeros({N, N})}
    {
    }

    /*!
     * @brief Perform Cholesky Factorization
     *
     * @param A Symmetric Matrix
     *
     * If $A$ is positive definite, then $p$ is zero.
     * If it is not, then $p$ is a positive integer,
     * such that $v = R^−1 e_p$ is a certificate vector
     * to make $v'*A[:p,:p]*v < 0$
     */
    void factorize(const Mat& A)
    {
        this->factor([&](unsigned i, unsigned j) { return A(i, j); });
    }

    /*!
     * @brief Perform Cholesky Factorization (Lazy evaluation)
     *
     * @tparam Fn
     * @param getA function to access the elements of A
     *
     * See also: factorize()
     */
    template <typename Fn>
    void factor(Fn getA)
    {
        auto& T = this->T;
        this->p = std::pair {0U, 0U};

        auto i = 0U;
        for (; i < this->n; ++i)
        {
            for (auto j = this->p.first; j <= i; ++j)
            {
                auto d = getA(i, j);
                for (auto k = this->p.first; k < j; ++k)
                {
                    d -= T(k, i) * T(j, k);
                }
                T(i, j) = d;
                if (i != j)
                {
                    T(j, i) = d / T(j, j);
                }
            }
            if constexpr (Allow_semidefinite)
            {
                if (T(i, i) < 0.)
                {
                    // this->stop = i + 1;
                    this->p.second = i + 1;
                    break;
                }
                if (T(i, i) == 0.)
                {
                    this->p.first = i + 1;
                }
            }
            else // strict positive definite
            {
                if (T(i, i) <= 0.)
                {
                    this->p.second = i + 1;
                    break;
                }
            }
        }
    }

    /*!
     * @brief Is $A$ symmetric positive definite (spd)
     *
     * @return true
     * @return false
     */
    auto is_spd() const -> bool
    {
        return this->p.second == 0;
    }

    /*!
     * @brief witness that certifies $A$ is not
     * symmetric positive definite (spd)
     *
     * @return auto
     */
    auto witness() -> double
    {
        if (this->is_spd())
        {
            throw std::runtime_error {"Implementation Error."};
        }
        auto [start, n] = this->p;
        auto m = n - 1; // assume stop >= 0
        this->v(m) = 1.;

        for (auto i = m; i > start; --i)
        {
            auto s = 0.;
            for (auto k = i; k <= m; ++k)
            {
                s += this->T(i - 1, k) * this->v(k);
            }
            this->v(i - 1) = -s;
        }
        return -this->T(m, m);
    }

    /*!
     * @brief Calculate v'*{A}(p,p)*v
     *
     * @param v
     * @param A
     * @return double
     */
    auto sym_quad(const Vec& A) const -> double
    {
        auto res = 0.;
        const auto& v = this->v;
        const auto& [start, stop] = this->p;
        for (auto i = start; i < stop; ++i)
        {
            auto s = 0.;
            for (auto j = i + 1; j < stop; ++j)
            {
                s += A(i, j) * v(j);
            }
            res += v(i) * (A(i, i) * v(i) + 2 * s);
        }
        return res;
    }

    auto sqrt() -> Mat
    {
        if (!this->is_spd())
        {
            throw std::runtime_error {"Implementation Error."};
        }

        // if (!this->sqrt_free) {
        //     return Mat{this->T};
        // }
        auto M = zeros({this->n, this->n});

        for (auto i = 0U; i < this->n; ++i)
        {
            M(i, i) = std::sqrt(this->T(i, i));
            for (auto j = i + 1; j < this->n; ++j)
            {
                M(i, j) = this->T(i, j) * M(i, i);
            }
        }

        return M;
    }
};
