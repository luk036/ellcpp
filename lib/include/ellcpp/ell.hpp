#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ELL_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ELL_HPP 1

#include <cmath>
#include <tuple>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

/* linux-2.6.38.8/include/linux/compiler.h */
#include <stdio.h>
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

/**
 * @brief Ellipsoid Search Space
 *
 * ell = { x | (x - xc)' * (\kappa Q)^-1 * (x - xc) <= 1 }
 */
class ell {
  public:
    using Arr = xt::xarray<double>;
    using params_t = std::tuple<double, double, double>;
    using return_t = std::tuple<int, params_t>;

  public:
    bool _use_parallel_cut = true;

  private:
    std::size_t _n;
    double _c1;
    double _kappa;
    Arr _xc;
    Arr _Q;

  public:
    /**
     * @brief Construct a new ell object
     *
     * @tparam T
     * @param val
     * @param x
     */
    template <typename T>
    ell(const T &val, const Arr &x)
        : _n{x.size()},                         // n
          _c1{double(_n * _n) / (_n * _n - 1)}, //
          _xc{x} {
        if constexpr (std::is_scalar<T>::value) { // C++17
            this->_Q = xt::eye(_n);
            this->_kappa = val;
        } else {
            this->_Q = xt::diag(val);
            this->_kappa = 1.;
        }
    }

    /**
     * @brief Construct a new ell object
     *
     * @param E
     */
    ell(const ell &E) = default;

    /**
     * @brief
     *
     * @return const Arr&
     */
    auto xc() const -> const Arr & { return _xc; }

    /**
     * @brief
     *
     * @return Arr&
     */
    auto xc() -> Arr & { return _xc; }

    /**
     * @brief Set the xc object
     *
     * @param xc
     */
    void set_xc(const Arr &xc) { _xc = xc; }

    /**
     * @brief
     *
     * @tparam T
     * @param g
     * @param beta
     * @return std::tuple<int, double>
     */
    template <typename T> //
    std::tuple<int, double> update(const Arr &g, const T &beta) {
        return this->update_core(g, beta);
    }

    /**
     * @brief Update ellipsoid core function using the cut
     *          g' * (x - xc) + beta <= 0
     *
     * @tparam T
     * @param g
     * @param beta
     * @return std::tuple<int, double>
     */
    template <typename T>
    std::tuple<int, double> update_core(const Arr &g, const T &beta) {
        auto Qg = Arr{xt::linalg::dot(_Q, g)};
        auto omega = xt::linalg::dot(g, Qg)();
        auto tsq = this->_kappa * omega;
        if (unlikely(tsq <= 0)) {
            return {4, 0.};
        }
        // auto tau = std::sqrt(_kappa * tsq);
        // auto alpha = beta / tau;
        auto [status, params] = this->calc_ll(beta, tsq);
        if (status != 0) {
            return {status, tsq};
        }
        const auto &[rho, sigma, delta] = params;
        this->_xc -= (rho / omega) * Qg;
        this->_Q -= (sigma / omega) * xt::linalg::outer(Qg, Qg);
        this->_kappa *= delta;
        if (unlikely(this->_kappa > 1e100 || this->_kappa < 1e-100)) {
            this->_Q *= this->_kappa;
            this->_kappa = 1.;
        }
        return {status, tsq}; // g++-7 is ok
    }

    /**
     * @brief parallel or deep cut
     *
     * @tparam T
     * @param beta
     * @param tsq
     * @return return_t
     */
    template <typename T> //
    return_t calc_ll(const T &beta, double tsq) {
        if constexpr (std::is_scalar<T>::value) { // C++17
            return this->calc_dc(beta, tsq);
        } else { // parallel cut
            if (unlikely(beta.shape()[0] < 2)) {
                return this->calc_dc(beta[0], tsq);
            }
            return this->calc_ll_core(beta[0], beta[1], tsq);
        }
    }

    /**
     * @brief
     *
     * @param b0
     * @param b1
     * @param tsq
     * @return return_t
     */
    return_t calc_ll_core(double b0, double b1, double tsq) const;

    /**
     * @brief Parallel Cut, one of them is central
     *
     * @param b1
     * @param t1
     * @param tsq
     * @return return_t
     */
    return_t calc_ll_cc(double b1, double t1, double tsq) const;

    /**
     * @brief Deep Cut
     *
     * @param b0
     * @param tsq
     * @return return_t
     */
    return_t calc_dc(double b0, double tsq) const;

    /**
     * @brief Central Cut
     *
     * @param tsq
     * @return return_t
     */
    return_t calc_cc(double tsq) const;
}; // } ell

/**
 * @brief Ellipsoid Method for special 1D case
 *
 */
class ell1d {
  public:
    using return_t = std::tuple<int, double>;

  private:
    double _r;
    double _xc;

  public:
    /**
     * @brief Construct a new ell1d object
     *
     * @param l
     * @param u
     */
    ell1d(double l, double u) //
        : _r{(u - l) / 2},    //
          _xc{l + _r} {}

    /**
     * @brief Construct a new ell1d object
     *
     * @param E
     */
    ell1d(const ell1d &E) = default;

    /**
     * @brief
     *
     * @return double
     */
    double xc() const { return _xc; }

    /**
     * @brief
     *
     * @param g
     * @param beta
     * @return return_t
     */
    return_t update(double g, double beta);
}; // } ell1d

#endif
