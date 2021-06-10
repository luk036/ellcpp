#include <cmath>
#include <ellcpp/cut_config.hpp>
#include <ellcpp/ell_assert.hpp>
#include <ellcpp/ell_stable.hpp>
// #include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

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
auto ell_stable::update(const std::tuple<Arr, T>& cut) -> std::tuple<CUTStatus, double>
{
    const auto& [g, beta] = cut;

    // calculate inv(L)*g: (n-1)*n/2 multiplications
    auto invLg = Arr {g}; // initially
    for (auto i = 1; i != this->_n; ++i)
    {
        for (auto j = 0; j != i; ++j)
        {
            this->_Q(i, j) = this->_Q(j, i) * invLg(j); 
            // keep for rank-one update
            invLg(i) -= this->_Q(i, j);
        }
    }

    // calculate inv(D)*inv(L)*g: n
    auto invDinvLg = Arr {invLg}; // initially
    for (auto i = 0; i != this->_n; ++i)
    {
        invDinvLg(i) *= this->_Q(i, i);
    }

    // calculate omega: n
    auto gQg = Arr {invDinvLg}; // initially
    auto omega = 0.;              // initially
    for (auto i = 0; i != this->_n; ++i)
    {
        gQg(i) *= invLg(i);
        omega += gQg(i);
    }

    this->_tsq = this->_kappa * omega;

    auto status = this->_update_cut(beta);
    if (status != CUTStatus::success)
    {
        return {status, this->_tsq};
    }

    // calculate Q*g = inv(L')*inv(D)*inv(L)*g : (n-1)*n/2
    auto Qg = Arr {invDinvLg}; // initially
    for (auto i = this->_n - 1; i > 0; --i)
    { // backward subsituition
        for (auto j = i; j != this->_n; ++j)
        {
            Qg(i - 1) -= this->_Q(i, j) * Qg(j); // ???
        }
    }

    // calculate xc: n
    this->_xc -= (this->_rho / omega) * Qg;

    // rank-one update: 3*n + (n-1)*n/2
    // const auto r = this->_sigma / omega;
    const auto mu = this->_sigma / (1. - this->_sigma);
    auto oldt = omega / mu; // initially
    const auto m = this->_n - 1;
    for (auto j = 0; j != m; ++j)
    {
        // p=sqrt(k)*vv(j);
        // const auto p = invLg(j);
        // const auto mup = mu * p;
        const auto t = oldt + gQg(j);
        // this->_Q(j, j) /= t; // update invD
        const auto beta2 = invDinvLg(j) / t;
        this->_Q(j, j) *= oldt / t; // update invD
        for (auto l = j + 1; l != this->_n; ++l)
        {
            // v(l) -= p * this->_Q(j, l);
            this->_Q(j, l) += beta2 * this->_Q(l, j);
        }
        oldt = t;
    }

    // const auto p = invLg(n1);
    // const auto mup = mu * p;
    const auto t = oldt + gQg(m);
    this->_Q(m, m) *= oldt / t; // update invD


    this->_kappa *= this->_delta;

    // if (this->no_defer_trick)
    // {
    //     this->_Q *= this->_kappa;
    //     this->_kappa = 1.;
    // }
    return {status, this->_tsq}; // g++-7 is ok
}

// Instantiation
template std::tuple<CUTStatus, double> ell_stable::update(
    const std::tuple<Arr, double>& cut);
template std::tuple<CUTStatus, double> ell_stable::update(
    const std::tuple<Arr, Arr>& cut);
