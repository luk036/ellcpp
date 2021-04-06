// -*- coding: utf-8 -*-
#pragma once

#include <ellcpp/utility.hpp>
#include <ellcpp/ell_assert.hpp> // ELL_UNLIKELY
#include <xtensor/xarray.hpp>

/*!
 * @brief Cholesky factorization for LMI
 *
 *  - LDL^T square-root-free version
 *  - Option allow semidefinite
 *  - A matrix A in R^{m x m} is positive definite iff v' A v > 0
 *      for all v in R^n.
 *  - O(p^2) per iteration, independent of N
 */
class chol_ext
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Vec = Arr;
    using Mat = Arr;
    using Rng = std::pair<size_t, size_t>;

  public:
    Rng p {0, 0}; //!< the rows where the process starts and stops
    Vec v;        //!< witness vector

  private:
    const size_t n; //!< dimension
    Mat T;          //!< temporary storage

  public:
    /*!
     * @brief Construct a new chol ext object
     *
     * @param[in] N dimension
     */
    explicit chol_ext(size_t N)
        : v {zeros({N})}
        , n {N}
        , T {zeros({N, N})}
    {
    }

    chol_ext(const chol_ext&) = delete;
    chol_ext& operator=(const chol_ext&) = delete;
    chol_ext(chol_ext&&) = default;

    /*!
     * @brief Perform Cholesky Factorization
     *
     * @param[in] A Symmetric Matrix
     *
     * If $A$ is positive definite, then $p$ is zero.
     * If it is not, then $p$ is a positive integer,
     * such that $v = R^-1 e_p$ is a certificate vector
     * to make $v'*A[:p,:p]*v < 0$
     */
    auto factorize(const Mat& A) -> bool
    {
        return this->factor([&](size_t i, size_t j) { return A(i, j); });
    }

    /*!
     * @brief Perform Cholesky Factorization (Lazy evaluation)
     *
     * @tparam Fn
     * @param[in] getA function to access the elements of A
     *
     * See also: factorize()
     */
    template <typename Callable, bool Allow_semidefinite = false>
    auto factor(Callable&& getA) -> bool
    {
        this->p = {0U, 0U};
        auto& [start, stop] = this->p;

        for (size_t i = 0U; i != this->n; ++i)
        {
            auto j = start;
            auto d = getA(i, j);
            while (j != i)
            {
                this->T(j, i) = d;
                this->T(i, j) = d / this->T(j, j); // note: T(j, i) here!
                d = getA(i, ++j);
                for (auto k = start; k != j; ++k)
                {
                    d -= this->T(i, k) * this->T(k, j);
                }
            }
            this->T(i, i) = d;

            if constexpr (Allow_semidefinite)
            {
                if (d < 0.)
                {
                    // this->stop = i + 1;
                    stop = i + 1;
                    break;
                }
                if (ELL_UNLIKELY(d == 0.))
                {
                    start = i + 1; 
                    // restart at i + 1, special as an LMI oracle
                }
            }
            else // not Allow_semidefinite
            {
                if (d <= 0.)
                {
                    stop = i + 1;
                    break;
                }
            }
        }

        return this->is_spd();
    }


    /*!
     * @brief Is $A$ symmetric positive definite (spd)
     *
     * @return true
     * @return false
     */
    auto is_spd() const noexcept -> bool
    {
        return this->p.second == 0;
    }

    /*!
     * @brief witness that certifies $A$ is not
     * symmetric positive definite (spd)
     *
     * @return auto
     */
    auto witness() -> double;

    /*!
     * @brief Calculate v'*{A}(p,p)*v
     *
     * @param[in] A
     * @return double
     */
    auto sym_quad(const Vec& A) const -> double;

    auto sqrt() -> Mat;
};
