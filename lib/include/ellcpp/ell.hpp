#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ELL_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ELL_HPP 1

#include <cmath>
#include <tuple>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

/* linux-2.6.38.8/include/linux/compiler.h */
#include <stdio.h>
#define likely(x)    __builtin_expect(!!(x), 1)
#define unlikely(x)  __builtin_expect(!!(x), 0)

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
    bool _use_parallel = true;

  private:
    std::size_t _n;
    double _c1;
    double _kappa;
    Arr _xc;
    Arr _Q;

  public:
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

    ell(const ell &E) = default;

    const Arr &xc() const { return _xc; }

    Arr &xc() { return _xc; }

    void set_xc(const Arr &xc) { _xc = xc; }

    template <typename T> //
    auto update(const Arr &g, const T &beta) {
        return this->update_core(g, beta);
    }

    /**
     * @brief Update ellipsoid core function using the cut
     *          g' * (x - xc) + beta <= 0
     *
     * @tparam T
     * @tparam V
     * @param g
     * @param beta
     * @return auto
     */
    template <typename T>
    std::tuple<int, double> update_core(const Arr &g, const T &beta) {
        auto Qg = xt::linalg::dot(_Q, g);
        auto omega = xt::linalg::dot(g, Qg)();
        auto tsq = this->_kappa * omega;
        if (unlikely(tsq <= 0.)) {
            return {4, 0.};
        }
        // auto tau = std::sqrt(_kappa * tsq);
        // auto alpha = beta / tau;
        auto [status, params] = this->calc_ll(beta, tsq);
        if (status != 0) {
            return {status, tsq};
        }
        auto [rho, sigma, delta] = params;
        this->_xc -= (rho / omega) * Qg;
        this->_Q -= (sigma / omega) * xt::linalg::outer(Qg, Qg);
        this->_kappa *= delta;
        return {status, tsq}; // g++-7 is ok
    }

    /* parallel or deep cut */
    template <typename T> //
    return_t calc_ll(const T &beta, double tsq) {
        if constexpr (std::is_scalar<T>::value) { // C++17
            return this->calc_dc(beta, tsq);
        } else { // parallel cut
            if (beta.shape()[0] < 2) {
                return this->calc_dc(beta[0], tsq);
            }
            return this->calc_ll_core(beta[0], beta[1], tsq);
        }
    }

    /**
     * @brief Core Parallel Cut
     *
     * @param a0
     * @param a1
     * @param n
     * @return auto
     */
    return_t calc_ll_core(double b0, double b1, double tsq) const;

    /**
     * @brief Parallel Cut, one of them is central
     *
     * @param a1
     * @param n
     * @return auto
     */
    return_t calc_ll_cc(double b1, double t1, double tsq) const;

    /**
     * @brief Deep Cut
     */
    return_t calc_dc(double b0, double tsq) const;

    /**
     * @brief Central Cut
     */
    return_t calc_cc(double tsq) const;
}; // } ell

class ell1d {
  public:
    using return_t = std::tuple<int, double>;

  private:
    double _r;
    double _xc;

  public:
    ell1d(double l, double u) //
        : _r{(u - l) / 2},    //
          _xc{l + _r} {}

    auto &xc() { return _xc; }

    ell1d(const ell1d &E) = default;

    return_t update(double g, double beta);
}; // } ell1d

#endif