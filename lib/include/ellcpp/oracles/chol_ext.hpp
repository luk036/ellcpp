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
    std::size_t _p;
    std::size_t _n;

  public:
    xt::xarray<double> _R{};

  public:
    /**
     * @brief Construct a new chol ext object
     *
     * @param n
     */
    explicit chol_ext(std::size_t n)
        : _p{0}, _n{n}, //
          _R{xt::zeros<double>({n, n})} {}

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
        _p = 0;

        for (auto i = 0u; i < _n; ++i) {
            for (auto j = 0u; j <= i; ++j) {
                d = getA(i, j);
                for (auto k = 0u; k < j; ++k) {
                    d -= _R(k, i) * _R(k, j);
                }
                if (i != j) {
                    _R(j, i) = d / _R(j, j);
                }
            }
            if (d <= 0) {
                _p = i + 1;
                _R(i, i) = std::sqrt(-d);
                break;
            }
            _R(i, i) = std::sqrt(d);
        }
    }

    /**
     * @brief
     *
     * @return true
     * @return false
     */
    bool is_spd() const { return _p == 0; }

    /**
     * @brief
     *
     * @return auto
     */
    auto witness() const {
        assert(!this->is_spd());
        Vec v = xt::zeros<double>({_p});
        using xt::placeholders::_;

        auto r = _R(_p - 1, _p - 1);
        auto ep = (r == 0) ? 0. : 1.;
        v[_p - 1] = (r == 0) ? 1. : 1. / r;

        for (int i = _p - 2; i >= 0; --i) {
            double s = 0.;
            for (auto k = i + 1; k < _p; ++k) {
                s += _R(i, k) * v[k];
            }
            v[i] = -s / _R(i, i);
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
        for (auto i = 0u; i < _p; ++i) {
            auto s = 0.;
            for (auto j = i + 1; j < _p; ++j) {
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

    for (; _p<_n; ++_p) {
        auto d = A(_p,_p);
        for (auto k=0u; k<_p; ++k ) {
            auto s = A(k,_p);
            for (auto i=0u; i<k; ++i ) {
                s -= _R(i,k)*_R(i,_p);
            }
            d -= sq( _R(k,_p) = s / _R(k,k) );
        }

        if (d < 0) {
            _R(_p,_p) = std::sqrt( -d );
            break;
        }
        _R(_p,_p) = sqrt( d );
    }
}

chol_ext::Vec chol_ext::witness() const
{
    chol_ext::Vec v(_p+1);
    v[_p] = 1.0 / _R(_p,_p);

    for (int i=_p-1; i>=0; --i) {
        auto s = 0.;
        for (int j=_p; j>i; --j) {
            s += _R(i,j) * v[j];
        }
        v[i] = -s / _R(i,i);
    }
    return v;
}

**/
#endif
