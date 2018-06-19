#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ELL_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ELL_HPP 1

#include <cmath>
#include <tuple>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

/**
 * @brief Ellipsoid Search Space
 *
 * ell = { x | (x - xc)' * P^-1 * (x - xc) <= 1 }
 */
class ell {
  public:
    using Arr = xt::xarray<double>;
    using params_t = std::tuple<double, double, double>;
    using return_t = std::tuple<int, params_t>;

  public:
    bool _use_parallel = true;

  private:
    std::size_t n;
    double _c1;
    double _kappa;
    Arr _xc;
    Arr _Q;

  public:
    template <typename T>
    ell(const T &val, const Arr &x)
        : n{x.size()} // n
        , _c1{n * n / (n * n - 1.)} // 
        , _xc{x} //
    {
        if constexpr (std::is_scalar<T>::value) { // C++17
            this->_Q = xt::eye(n);
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
    auto update_core(const Arr &g, const T &beta) {
        auto Qg = xt::linalg::dot(_Q, g);
        auto tsq = xt::linalg::dot(g, Qg)();
        if (tsq <= 0.) {
            return std::tuple{4, 0.};
        }
        auto tau = std::sqrt(_kappa * tsq);
        auto alpha = beta / tau;
        auto [status, params] = this->calc_ll(alpha);
        if (status != 0) {
            return std::tuple{status, tau}; // g++-7 is ok with this
        }
        auto [rho, sigma, delta] = params;
        this->_xc -= (_kappa * rho / tau) * Qg;
        this->_Q -= (sigma / tsq) * xt::linalg::outer(Qg, Qg);
        this->_kappa *= delta;
        return std::tuple{status, tau}; // g++-7 is ok
    }

    /**
     * @brief Central Cut
     */
    params_t calc_cc(std::size_t n) const;

    /**
     * @brief Deep Cut
     */
    return_t calc_dc(double alpha) const;

    /**
     * @brief Parallel Cut, one of them is central
     *
     * @param a1
     * @param n
     * @return auto
     */
    params_t calc_ll_cc(double a1, std::size_t n) const;

    /**
     * @brief General Parallel Cut
     *
     * @param a0
     * @param a1
     * @param n
     * @return auto
     */
    params_t calc_ll_general(double a0, double a1, std::size_t n) const;

    /* parallel or deep cut */
    template <typename T> auto calc_ll(const T &alpha) {
        if constexpr (std::is_scalar<T>::value) { // C++17
            return this->calc_dc(alpha);
        } else { // parallel cut
            auto a0 = alpha[0];
            if (alpha.shape()[0] < 2) {
                return this->calc_dc(a0);
            }
            auto a1 = alpha[1];
            if (a1 >= 1. || !_use_parallel) {
                return this->calc_dc(a0);
            }
            auto n = _xc.size();
            auto status = 0;
            auto params = std::tuple{0., 0., 0.};
            
            if (a0 > a1) {
                status = 1; // no sol'n
            } else if (n*a0*a1 < -1.) {
                status = 3; // no effect
            } else if (a0 == 0.) {
                params = this->calc_ll_cc(a1, n);
            } else {
                params = this->calc_ll_general(a0, a1, n);
            }
            return std::tuple{status, params};
        }
    }

    template <typename T> auto update(const Arr &g, const T &beta) {
        return this->update_core(g, beta);
    }

}; // } ell

class ell1d {
  public:
    using return_t = std::tuple<int, double>;

  private:
    double _r;
    double _xc;

  public:
    ell1d(double l, double u) : _r{(u - l) / 2}, _xc{l + _r} {}

    auto &xc() { return _xc; }

    ell1d(const ell1d &E) = default;

    return_t update(double g, double beta);
}; // } ell1d

#endif