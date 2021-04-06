#include <cmath>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/ell_assert.hpp>
#include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;


/*!
 * @brief
 *
 * @param[in] b0
 * @param[in] b1
 * @return int
 */
CUTStatus ell::_calc_ll_core(const double& b0, const double& b1)
{
    //const auto b1sq = b1 * b1;
    const auto b1sqn = b1 * (b1 / this->_tsq);
    const auto t1n = 1. - b1sqn;
    if (t1n < 0. || !this->use_parallel_cut)
    {
        return this->_calc_dc(b0);
    }

    const auto bdiff = b1 - b0;
    if (bdiff < 0.)
    {
        return CUTStatus::nosoln; // no sol'n
    }

    if (b0 == 0.) // central cut
    {
        this->_calc_ll_cc(b1, b1sqn);
        return CUTStatus::success;
    }

    const auto b0b1 = b0 * b1;
    const auto& n = this->_n;
    if (ELL_UNLIKELY(n * b0b1 < -this->_tsq))
    {
        return CUTStatus::noeffect; // no effect
    }

    // const auto t0 = this->_tsq - b0 * b0;
    const auto t0n = 1. - b0 * (b0 / this->_tsq);
    // const auto t1 = this->_tsq - b1sq;
    const auto bav = (b0 + b1) / 2;
    const auto bavn = bav / this->_tsq;
    const auto tempn = n * bavn * (b1 - b0);
    const auto xi = std::sqrt(t0n * t1n + tempn * tempn);
    this->_sigma = (n + (1. - b0b1 / this->_tsq - xi) / (2 * bavn * bav)) / this->_nPlus1;
    this->_rho = this->_sigma * bav;
    this->_delta = this->_c1 * ((t0n + t1n) / 2. + xi / n);
    return CUTStatus::success;
}

/*!
 * @brief
 *
 * @param[in] b1
 * @param[in] b1sq
 * @return void
 */
void ell::_calc_ll_cc(const double& b1, const double& b1sqn)
{
    const auto temp = this->_n * b1sqn / 2.;
    const auto xi = std::sqrt((1. - b1sqn) + temp * temp);
    this->_sigma = (this->_n + 2. * (1. - xi) / b1sqn) / this->_nPlus1;
    this->_rho = this->_sigma * b1 / 2;
    this->_delta = this->_c1 * (1. - b1sqn / 2.  + xi / this->_n);
}

/*!
 * @brief Deep Cut
 *
 * @param[in] beta
 * @return int
 */
CUTStatus ell::_calc_dc(const double& beta)
{
    const auto tau = std::sqrt(this->_tsq);

    const auto bdiff = tau - beta;
    if (bdiff < 0)
    {
        return CUTStatus::nosoln; // no sol'n
    }

    if (beta == 0.)
    {
        this->_calc_cc(tau);
        return CUTStatus::success;
    }

    const auto gamma = tau + this->_n * beta;
    if (ELL_UNLIKELY(gamma < 0))
    {
        return CUTStatus::noeffect; // no effect
    }

    this->_mu = (bdiff / gamma) * this->_nMinus1 / 2.;
    this->_rho = gamma / (this->_n + 1.);
    this->_sigma = 2 * this->_rho / (tau + beta);
    this->_delta = this->_c1 * (this->_tsq - beta * beta) / this->_tsq;
    return CUTStatus::success;
}

/*!
 * @brief Central Cut
 *
 * @param[in] tau
 * @return int
 */
void ell::_calc_cc(const double& tau)
{
    this->_mu = this->_nMinus1 / 2.;
    this->_sigma = 2. / this->_nPlus1;
    this->_rho = tau / this->_nPlus1;
    this->_delta = this->_c1;
}


/*!
 * @brief Update ellipsoid core function using the cut
 *
 *        g' * (x - xc) + beta <= 0
 *
 * @tparam T
 * @param[in] cut
 * @return std::tuple<int, double>
 */
template <typename T>
std::tuple<CUTStatus, double> ell::update(const std::tuple<Arr, T>& cut)
{
    // const auto& [g, beta] = cut;
    const auto& beta = std::get<1>(cut);

    const auto& g = std::get<0>(cut);
    // n^2
    const auto Qg = Arr {xt::linalg::dot(this->_Q, g)}; // n^2
    const auto omega = xt::linalg::dot(g, Qg)();        // n
    this->_tsq = this->_kappa * omega;

    auto status = this->_update_cut(beta);
    if (status != CUTStatus::success)
    {
        return {status, this->_tsq};
    }

    this->_xc -= (this->_rho / omega) * Qg; // n
    // n*(n+1)/2 + n
    // this->_Q -= (this->_sigma / omega) * xt::linalg::outer(Qg, Qg);
    const auto r = this->_sigma / omega;
    for (auto i = 0U; i < this->_n; ++i)
    {
        const auto rQg = r * Qg(i);
        for (auto j = 0U; j < i; ++j)
        {
            this->_Q(i, j) -= rQg * Qg(j);
            this->_Q(j, i) = this->_Q(i, j);
        }
        this->_Q(i, i) -= rQg * Qg(i);
    }

    this->_kappa *= this->_delta;

    if (this->no_defer_trick)
    {
        this->_Q *= this->_kappa;
        this->_kappa = 1.;
    }
    return {status, this->_tsq}; // g++-7 is ok
}

// Instantiation
template std::tuple<CUTStatus, double> ell::update(
    const std::tuple<Arr, double>& cut);
template std::tuple<CUTStatus, double> ell::update(
    const std::tuple<Arr, Arr>& cut);
