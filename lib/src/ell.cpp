#include <cmath>
#include <ellcpp/ell.hpp>
// #include <tuple>
#include <xtensor-blas/xlinalg.hpp>

/* linux-2.6.38.8/include/linux/compiler.h */
#include <cstdio>
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

using Arr = xt::xarray<double, xt::layout_type::row_major>;

/*!
 * @brief
 *
 * @param b0
 * @param b1
 * @param tsq
 * @return ell::return_t
 */
ell::return_t ell::__calc_ll_core(double b0, double b1, double tsq) const
{
    auto b1sq = b1 * b1;
    if (b1sq > tsq || !this->_use_parallel_cut)
    {
        return this->__calc_dc(b0, tsq);
    }

    auto params = std::tuple {0., 0., 0.};

    if (unlikely(b1 < b0))
    {
        return {1, std::move(params)}; // no sol'n
    }

    if (b0 == 0.)
    {
        return this->__calc_ll_cc(b1, b1sq, tsq);
    }

    auto n = this->_n;
    auto b0b1 = b0 * b1;
    if (unlikely(n * b0b1 < -tsq))
    {
        return {3, std::move(params)}; // no effect
    }

    auto t0 = tsq - b0 * b0;
    auto t1 = tsq - b1sq;
    auto bav = (b0 + b1) / 2;
    auto temp = n * bav * (b1 - b0);
    auto xi = std::sqrt(t0 * t1 + temp * temp);
    auto sigma = (n + (tsq - b0b1 - xi) / (2 * bav * bav)) / (n + 1.);
    auto rho = sigma * bav;
    auto delta = this->_c1 * ((t0 + t1) / 2 + xi / n) / tsq;
    params = std::tuple {rho, sigma, delta};
    return {0, std::move(params)};
}

/*!
 * @brief
 *
 * @param b1
 * @param b1sq
 * @param tsq
 * @return ell::return_t
 */
ell::return_t ell::__calc_ll_cc(double b1, double b1sq, double tsq) const
{
    auto n = this->_n;
    auto temp = n * b1sq / 2;
    auto xi = std::sqrt(tsq * (tsq - b1sq) + temp * temp);
    auto sigma = (n + 2 * (tsq - xi) / b1sq) / (n + 1.);
    auto rho = sigma * b1 / 2;
    auto delta = this->_c1 * (tsq - b1sq / 2 + xi / n) / tsq;
    auto params = ell::params_t {rho, sigma, delta};
    return {0, std::move(params)};
}

/*!
 * @brief Deep Cut
 *
 * @param beta
 * @param tsq
 * @return ell::return_t
 */
ell::return_t ell::__calc_dc(double beta, double tsq) const
{
    auto params = std::tuple {0., 0., 0.};
    auto tau = std::sqrt(tsq);

    if (beta > tau)
    {
        return {1, std::move(params)}; // no sol'n
    }

    if (beta == 0.)
    {
        return this->__calc_cc(tsq);
    }

    auto n = this->_n;
    auto gamma = tau + n * beta;
    if (unlikely(gamma < 0))
    {
        return {3, std::move(params)}; // no effect
    }

    auto rho = gamma / (n + 1.);
    auto sigma = 2 * rho / (tau + beta);
    auto delta = this->_c1 * (tsq - beta * beta) / tsq;
    auto ret = ell::params_t {rho, sigma, delta};
    return {0, std::move(ret)};
}

/*!
 * @brief Central Cut
 *
 * @param tsq
 * @return ell::return_t
 */
ell::return_t ell::__calc_cc(double tsq) const
{
    auto np1 = this->_n + 1;
    auto sigma = 2. / np1;
    auto rho = std::sqrt(tsq) / np1;
    auto delta = this->_c1;
    auto params = ell::params_t {rho, sigma, delta};
    return {0, std::move(params)};
}

/*!
 * @brief
 *
 * @param g
 * @param beta
 * @return ell1d::return_t
 */
ell1d::return_t ell1d::update(double g, double beta)
{
    auto tau = std::abs(this->_r * g);
    auto tsq = tau * tau;

    if (beta == 0.)
    {
        this->_r /= 2;
        this->_xc += g > 0. ? -this->_r : this->_r;
        return {0, tsq};
    }
    if (beta > tau)
    {
        return {1, tsq}; // no sol'n
    }
    if (unlikely(beta < -tau))
    {
        return {3, tsq}; // no effect
    }

    auto bound = this->_xc - beta / g;
    auto u = g > 0. ? bound : this->_xc + this->_r;
    auto l = g > 0. ? this->_xc - this->_r : bound;

    this->_r = (u - l) / 2;
    this->_xc = l + this->_r;
    return {0, tsq};
}

/*!
 * @brief Update ellipsoid core function using the cut
 *          g' * (x - xc) + beta <= 0
 *
 * @tparam T
 * @param g
 * @param beta
 * @return std::tuple<int, double>
 */
template <typename T>
std::tuple<int, double> ell::update(const Arr& g, const T& beta)
{
    auto Qg = Arr {xt::linalg::dot(_Q, g)};
    auto omega = xt::linalg::dot(g, Qg)();
    auto tsq = this->_kappa * omega;
    auto status = 0;
    auto params = std::tuple {0., 0., 0.};

    if constexpr (std::is_scalar<T>::value)
    { // C++17
        std::tie(status, params) = this->__calc_dc(beta, tsq);
    }
    else
    { // parallel cut
        if (unlikely(beta.shape()[0] < 2))
        {
            std::tie(status, params) = this->__calc_dc(beta[0], tsq);
        }
        else
        {
            std::tie(status, params) =
                this->__calc_ll_core(beta[0], beta[1], tsq);
        }
    }

    if (status != 0)
    {
        return {status, tsq};
    }
    auto& [rho, sigma, delta] = params;

    this->_xc -= (rho / omega) * Qg;
    this->_Q -= (sigma / omega) * xt::linalg::outer(Qg, Qg);
    this->_kappa *= delta;

    // if (unlikely(this->_kappa > 1e100 || this->_kappa < 1e-100)) {
    //     this->_Q *= this->_kappa;
    //     this->_kappa = 1.;
    // }
    return {status, tsq}; // g++-7 is ok
}

// Instantiation
template std::tuple<int, double> ell::update(const Arr& g, const double& beta);
template std::tuple<int, double> ell::update(const Arr& g, const Arr& beta);
