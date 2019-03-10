#ifndef CHOL_EXT_HPP
#define CHOL_EXT_HPP 1

// #include <boost/numeric/ublas/symmetric.hpp>
// #include <boost/numeric/ublas/triangular.hpp>
// #include <cassert>

// namespace bnu = boost::numeric::ublas;

#include <cassert>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

#include <iostream>

/**
 * @brief Cholesky factorization
 *
 */
class chol_ext {
    // using Mat = bnu::symmetric_matrix<double, bnu::upper>;
    // using UTMat = bnu::triangular_matrix<double, bnu::upper>;
    // using Vec = bnu::vector<double>;
    using Vec = xt::xarray<double>;
    using Mat = xt::xarray<double>;
    using shape_type = Vec::shape_type;

  private:
    std::size_t p;
    std::size_t n;

  public:
    xt::xarray<double> R{};

  public:
    /**
     * @brief Construct a new chol ext object
     *
     * @param n
     */
    explicit chol_ext(std::size_t n)
        : p{0}, n{n}, //
          R{xt::zeros<double>({n, n})} {}

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
     * @brief
     *
     * @tparam Fn
     * @param getA
     *
     * If $A$ is positive definite, then $p$ is zero.
     * If it is not, then $p$ is a positive integer,
     * such that $v = R^{-1} e_p$ is a certificate vector
     * to make $v'*A[:p,:p]*v < 0$
     */
    template <typename Fn> void factor(Fn getA) {
        double d;
        this->p = 0;
        auto &R = this->R;

        for (auto i = 0U; i < this->n; ++i) {
            for (auto j = 0U; j <= i; ++j) {
                d = getA(i, j);
                for (auto k = 0U; k < j; ++k) {
                    d -= R(k, i) * R(k, j);
                }
                if (i != j) {
                    R(j, i) = d / R(j, j);
                }
            }
            if (d <= 0) {
                this->p = i + 1;
                R(i, i) = std::sqrt(-d);
                break;
            }
            R(i, i) = std::sqrt(d);
        }
    }

    /**
     * @brief
     *
     * @return true
     * @return false
     */
    bool is_spd() const { return this->p == 0; }

    /**
     * @brief
     *
     * @return auto
     */
    auto witness() const {
        assert(!this->is_spd());
        auto &p = this->p;
        Vec v = xt::zeros<double>({p});
        using xt::placeholders::_;

        auto r = this->R(p - 1, p - 1);
        auto ep = (r == 0) ? 0. : 1.;
        v[p - 1] = (r == 0) ? 1. : 1. / r;

        for (int i = p - 2; i >= 0; --i) {
            double s = 0.;
            for (auto k = i + 1; k < p; ++k) {
                s += this->R(i, k) * v[k];
            }
            v[i] = -s / this->R(i, i);
        }
        return std::tuple{std::move(v), ep};
    }

    /**
     * @brief
     *
     * @param v
     * @param A
     * @return double
     */
    double sym_quad(const xt::xarray<double> &v, const xt::xarray<double> &A) {
        auto res = 0.;
        for (auto i = 0U; i < this->p; ++i) {
            auto s = 0.;
            for (auto j = i + 1; j < this->p; ++j) {
                s += A(i, j) * v(j);
            }
            res += v(i) * (A(i, i) * v(i) + 2 * s);
            // res += v(i) * s;
        }
        return res;
    }
};

/**
 * Constructs a upper triangular matrix R, such that R'*R= A.
 * If A is not symmetric positive-definite (SPD), only a partial
 * factorization is performed. If isspd() evalutate true then
 * the factorizaiton was successful.

chol_ext::chol_ext(const chol_ext::Mat& A )
    : _p{0}, _n{A.size1()}, _R(_n, _n) {
    auto sq = [](auto a) { return a*a; }; // square

    for (; _p<this->n; ++_p) {
        auto d = A(_p,_p);
        for (auto k=0U; k<_p; ++k ) {
            auto s = A(k,_p);
            for (auto i=0U; i<k; ++i ) {
                s -= this->R(i,k)*this->R(i,_p);
            }
            d -= sq( this->R(k,_p) = s / this->R(k,k) );
        }

        if (d < 0) {
            this->R(_p,_p) = std::sqrt( -d );
            break;
        }
        this->R(_p,_p) = sqrt( d );
    }
}

chol_ext::Vec chol_ext::witness() const
{
    chol_ext::Vec v(_p+1);
    v[_p] = 1.0 / this->R(_p,_p);

    for (int i=_p-1; i>=0; --i) {
        auto s = 0.;
        for (int j=_p; j>i; --j) {
            s += this->R(i,j) * v[j];
        }
        v[i] = -s / this->R(i,i);
    }
    return v;
}

**/
#endif
