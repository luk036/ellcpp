// -*- coding: utf-8 -*-
#pragma once

#include <stdexcept>
#include <cassert>
#include <xtensor/xarray.hpp>

/*!
 * @brief Cholesky factorization
 */
template<bool Allow_semidefinite = false> //
class chol_ext
{
    // using Mat = bnu::symmetric_matrix<double, bnu::upper>;
    // using UTMat = bnu::triangular_matrix<double, bnu::upper>;
    // using Vec = bnu::vector<double>;
    using Vec = xt::xarray<double, xt::layout_type::row_major>;
    using Mat = xt::xarray<double, xt::layout_type::row_major>;

public:
    std::size_t start;
    std::size_t stop;
    Vec         v;

private:
    std::size_t n;
    Mat         T;

public:
    /*!
     * @brief Construct a new chol ext object
     *
     * @param n
     */
    explicit chol_ext(std::size_t N) : v{xt::zeros<double>({N})}, n{N}, T{xt::zeros<double>({N, N})}
    {
    }

    /*!
     * @brief
     *
     * @param A
     *
     * If $A$ is positive definite, then $p$ is zero.
     * If it is not, then $p$ is a positive integer,
     * such that $v = R^{-1} e_p$ is a certificate vector
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
     * @param getA
     *
     * See also: factorize()
     */
    template<typename Fn>
    void factor(Fn getA)
    {
        auto& T     = this->T;
        auto  i     = 0U;
        this->start = 0U;
        this->stop  = 0U;

        for (; i < this->n; ++i)
        {
            for (auto j = this->start; j <= i; ++j)
            {
                auto d = getA(i, j);
                for (auto k = this->start; k < j; ++k)
                {
                    d -= T(k, i) * T(j, k);
                }
                T(i, j) = d;
                if (i != j) { T(j, i) = d / T(j, j); }
            }
            if constexpr (Allow_semidefinite)
            {
                if (T(i, i) < 0.)
                {
                    this->stop = i + 1;
                    break;
                }
                if (T(i, i) == 0.) { this->start = i + 1; }
            }
            else
            {
                if (T(i, i) <= 0.)
                {
                    this->stop = i + 1;
                    break;
                }
            }
            // if (T(i, i) > 0.) {
            //     continue;
            // }
            // if (T(i, i) < 0. || !this->allow_semidefinite) {
            //     break;
            // }
            // this->start = i+1;
        }
        // this->p = i;
    }

    /*!
     * @brief Is $A$ symmetric positive definite (spd)
     *
     * @return true
     * @return false
     */
    auto is_spd() const -> bool { return this->stop == 0; }

    /*!
     * @brief witness that certifies $A$ is not
     * symmetric positive definite (spd)
     *
     * @return auto
     */
    auto witness() -> double
    {
        if (this->is_spd()) { throw std::runtime_error{"Implementation Error."}; }
        // auto &p = this->p;
        auto&  stop = this->stop;
        size_t p    = stop - 1; // assume stop >= 1
        this->v(p)  = 1.;

        for (auto i = p; i > this->start; --i)
        {
            auto s = 0.;
            for (auto k = i; k <= p; ++k)
            {
                s += this->T(i - 1, k) * this->v(k);
            }
            this->v(i - 1) = -s;
        }

        return -this->T(p, p);
    }

    /*!
     * @brief
     *
     * @param v
     * @param A
     * @return double
     */
    double sym_quad(const Vec& A) const
    {
        auto  res = 0.;
        auto& v   = this->v;
        for (auto i = this->start; i < this->stop; ++i)
        {
            auto s = 0.;
            for (auto j = i + 1; j < this->stop; ++j)
            {
                s += A(i, j) * v(j);
            }
            res += v(i) * (A(i, i) * v(i) + 2 * s);
        }
        return res;
    }

    auto sqrt() -> Mat
    {
        if (!this->is_spd()) { throw std::runtime_error{"Implementation Error."}; }

        // if (!this->sqrt_free) {
        //     return Mat{this->T};
        // }
        auto M = Mat{xt::zeros<double>({this->n, this->n})};

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
