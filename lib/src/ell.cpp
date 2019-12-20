#include <cmath>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

/*!
 * @brief
 *
 * @param b0
 * @param b1
 * @return int
 */
CUTStatus ell::_calc_ll_core(const double& b0, const double& b1)
{
    const auto b1sq = b1 * b1;
    if (b1sq > this->_tsq or not this->_use_parallel_cut)
    {
        return this->_calc_dc(b0);
    }

    [[unlikely]] if (b1 < b0)
    {
        return CUTStatus::nosoln; // no sol'n
    }

    if (b0 == 0.)
    {
        this->_calc_ll_cc(b1, b1sq);
        return CUTStatus::success;
    }

    const auto b0b1 = b0 * b1;
    const auto& n = this->_n;
    [[unlikely]] if (n * b0b1 < -this->_tsq)
    {
        return CUTStatus::noeffect; // no effect
    }

    const auto t0 = this->_tsq - b0 * b0;
    const auto t1 = this->_tsq - b1sq;
    const auto bav = (b0 + b1) / 2;
    const auto temp = n * bav * (b1 - b0);
    const auto xi = std::sqrt(t0 * t1 + temp * temp);
    this->_sigma = (n + (this->_tsq - b0b1 - xi) / (2 * bav * bav)) / (n + 1.);
    this->_rho = this->_sigma * bav;
    this->_delta = this->_c1 * ((t0 + t1) / 2 + xi / n) / this->_tsq;
    return CUTStatus::success;
}

/*!
 * @brief
 *
 * @param b1
 * @param b1sq
 * @return void
 */
void ell::_calc_ll_cc(const double& b1, const double& b1sq)
{
    const auto& n = this->_n;
    const auto temp = n * b1sq / 2;
    const auto xi = std::sqrt(this->_tsq * (this->_tsq - b1sq) + temp * temp);
    this->_sigma = (n + 2 * (this->_tsq - xi) / b1sq) / (n + 1);
    this->_rho = this->_sigma * b1 / 2;
    this->_delta = this->_c1 * (this->_tsq - b1sq / 2 + xi / n) / this->_tsq;
}

/*!
 * @brief Deep Cut
 *
 * @param beta
 * @return int
 */
CUTStatus ell::_calc_dc(const double& beta)
{
    const auto tau = std::sqrt(this->_tsq);

    if (beta > tau)
    {
        return CUTStatus::nosoln; // no sol'n
    }

    if (beta == 0.)
    {
        this->_calc_cc(tau);
        return CUTStatus::success;
    }

    const auto gamma = tau + this->_n * beta;
    [[unlikely]] if (gamma < 0)
    {
        return CUTStatus::noeffect; // no effect
    }

    this->_rho = gamma / (this->_n + 1.);
    this->_sigma = 2 * this->_rho / (tau + beta);
    this->_delta = this->_c1 * (this->_tsq - beta * beta) / this->_tsq;
    return CUTStatus::success;
}

/*!
 * @brief Central Cut
 *
 * @param tau
 * @return int
 */
void ell::_calc_cc(const double& tau)
{
    const auto np1 = this->_n + 1;
    this->_sigma = 2. / np1;
    this->_rho = tau / np1;
    this->_delta = this->_c1;
}

/*!
 * @brief Update ellipsoid core function using the cut
 *
 *        g' * (x - xc) + beta <= 0
 *
 * @tparam T
 * @param g
 * @param beta
 * @return std::tuple<int, double>
 */
template <typename T>
std::tuple<CUTStatus, double> ell::update(const std::tuple<Arr, T>& cut)
{
    const auto& [g, beta] = cut;
    const auto Qg = Arr {xt::linalg::dot(_Q, g)};
    const auto omega = xt::linalg::dot(g, Qg)();
    this->_tsq = this->_kappa * omega;
    auto status = CUTStatus::success;

    if constexpr (std::is_scalar_v<T>)
    { // C++17
        status = this->_calc_dc(beta);
    }
    else
    { // parallel cut
        [[unlikely]] if (beta.shape()[0] < 2)
        {
            status = this->_calc_dc(beta[0]);
        }
        else
        {
            status = this->_calc_ll_core(beta[0], beta[1]);
        }
    }

    if (status != CUTStatus::success)
    {
        return {status, this->_tsq};
    }

    this->_xc -= (this->_rho / omega) * Qg;
    this->_Q -= (this->_sigma / omega) * xt::linalg::outer(Qg, Qg);
    this->_kappa *= this->_delta;

    if (this->_no_defer_trick)
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
