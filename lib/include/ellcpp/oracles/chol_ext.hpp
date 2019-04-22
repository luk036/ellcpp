#ifndef CHOL_EXT_HPP
#define CHOL_EXT_HPP 1

#include <cassert>
#include <xtensor/xarray.hpp>
#include <stdexcept>

/**
 * @brief Cholesky factorization
 */
class chol_ext {
    // using Mat = bnu::symmetric_matrix<double, bnu::upper>;
    // using UTMat = bnu::triangular_matrix<double, bnu::upper>;
    // using Vec = bnu::vector<double>;
    using Vec = xt::xarray<double>;
    using Mat = xt::xarray<double>;

  private:

    bool sqrt_free = true;
    std::size_t p;
    std::size_t n;
    Mat T;

  public:

    /**
     * @brief Construct a new chol ext object
     *
     * @param n
     */
    explicit chol_ext(std::size_t N)
        : n{N}, T{xt::zeros<double>({N, N})} {}

    /**
     * @brief
     *
     * @param A
     *
     * If $A$ is positive definite, then $p$ is zero.
     * If it is not, then $p$ is a positive integer,
     * such that $v = R^{-1} e_p$ is a certificate vector
     * to make $v'*A[:p,:p]*v < 0$
     */
    void factorize(const Mat &A) {
        this->factor([&](unsigned i, unsigned j) { return A(i, j); });
    }

    /**
     * @brief Perform Cholesky Factorization (Lazy evaluation)
     *
     * @tparam Fn
     * @param getA
     *
     * See also: factorize()
     */
    template <typename Fn> void factor(Fn getA) {
        auto& T = this->T;
        auto i = 0U;
        for (; i < this->n; ++i) {
            for (auto j = 0U; j <= i; ++j) {
                auto d = getA(i, j);
                for (auto k = 0U; k < j; ++k) {
                    d -= T(k, i) * T(j, k);
                }
                T(i, j) = d;
                if (i != j) {
                    T(j, i) = d / T(j, j);
                }
            }
            if (T(i, i) <= 0.) {
                break;
            }
        }
        this->p = i;
    }

    /**
     * @brief Perform Cholesky Factorization (Lazy evaluation)
     *
     * @tparam Fn
     * @param getA
     *
     * See also: factorize()
     */
    // template <typename Fn> void factor3(Fn getA) {
    //     this->sqrt_free = false;
    //     this->p = 0;

    //     for (auto i = 0U; i < this->n; ++i) {
    //         double d;
    
    //         for (auto j = 0U; j <= i; ++j) {
    //             d = getA(i, j);
    //             for (auto k = 0U; k < j; ++k) {
    //                 d -= this->T(k, i) * this->T(k, j);
    //             }
    //             if (i != j) {
    //                 this->T(j, i) = d / this->T(j, j);
    //             }
    //         }
    //         if (d <= 0.) {
    //             this->p = i + 1;
    //             this->T(i, i) = std::sqrt(-d);
    //             break;
    //         } else {
    //             this->T(i, i) = std::sqrt(d);
    //         }
    //     }
    // }

    /**
     * @brief Is $A$ symmetric positive definite (spd)
     *
     * @return true
     * @return false
     */
    auto is_spd() const -> bool { return this->p == this->n; }

    /**
     * @brief witness that certifies $A$ is not
     * symmetric positive definite (spd)
     *
     * @return auto
     */
    auto witness() const {
        if (this->is_spd()) {
            throw std::runtime_error{"Implementation Error."};
        }
        auto &p = this->p;
        auto n = p + 1;
        auto v = Vec{xt::zeros<double>({n})};
        // auto r = this->T(p - 1, p - 1);
        // auto ep = (r == 0.) ? 0. : 1.;
        // v[p - 1] = (r == 0.) ? 1. : 1. / r;
        v[p] = 1.;

        for (int i = p; i > 0; --i) {
            auto s = 0.;
            for (auto k = i; k < n; ++k) {
                s += this->T(i-1, k) * v[k];
            }
            v[i-1] = -s;
        }

        return std::tuple{std::move(v), -this->T(p, p)};
    }

    /**
     * @brief witness that certifies $A$ is not
     * symmetric positive definite (spd)
     *
     * @return auto
     */
    // auto witness3() const {
    //     if (this->sqrt_free) {
    //         throw std::runtime_error{"Implementation Error."};
    //     }
    //     if (this->is_spd()) {
    //         throw std::runtime_error{"Implementation Error."};
    //     }
    //     auto &p = this->p;
    //     auto v = Vec{xt::zeros<double>({p})};
    //     // auto r = this->T(p - 1, p - 1);
    //     // auto ep = (r == 0.) ? 0. : 1.;
    //     // v[p - 1] = (r == 0.) ? 1. : 1. / r;
    //     v[p - 1] = 1.;

    //     for (int i = p - 2; i >= 0; --i) {
    //         double s = 0.;
    //         for (auto k = i + 1; k < p; ++k) {
    //             s += this->T(i, k) * v[k];
    //         }
    //         v[i] = -s / this->T(i, i);
    //     }

    //     auto ep = this->T(p - 1, p - 1);
    //     return std::tuple{std::move(v), ep * ep};
    // }

    auto sqrt() -> Mat {
        if (!this->is_spd()) {
            throw std::runtime_error{"Implementation Error."};
        }

        if (!this->sqrt_free) {
            return Mat{this->T};
        }

        auto n = this->n;
        auto M = Mat{xt::zeros<double>({n, n})};

        for (auto i = 0U; i < n; ++i) {
            M(i, i) = std::sqrt(this->T(i, i));
            for (auto j = i+1; j < n; ++j) {
                M(i, j) = this->T(i, j) * M(i, i); 
            }
        }

        return std::move(M);
    }

    /**
     * @brief
     *
     * @param v
     * @param A
     * @return double
     */
    double sym_quad(const Vec &v, const Vec &A) {
        auto res = 0.;
        auto n = this->p + 1;
        for (auto i = 0U; i < n; ++i) {
            auto s = 0.;
            for (auto j = i + 1; j < n; ++j) {
                s += A(i, j) * v(j);
            }
            res += v(i) * (A(i, i) * v(i) + 2 * s);
        }
        return res;
    }
};

#endif
