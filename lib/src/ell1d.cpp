#include <cmath>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell1d.hpp>
#include <ellcpp/half_nonnegative.hpp>

// #include <tuple>

/*!
 * @brief
 *
 * @param[in] cut
 * @return ell1d::return_t
 */
ell1d::return_t ell1d::update(const std::tuple<double, double>& cut)
{
    const auto& [g, beta] = cut;

    const auto tau = std::abs(this->_r * g);
    const auto tsq = tau * tau;

    if (beta == 0.)
        [[unlikely]]
        {
            this->_r /= 2;
            this->_xc += g > 0. ? -this->_r : this->_r;
            return {CUTStatus::success, tsq};
        }
    if (beta > tau)
        [[unlikely]]
        {
            return {CUTStatus::nosoln, tsq}; // no sol'n
        }
    if (beta < -tau)
        [[unlikely]]
        {
            return {CUTStatus::noeffect, tsq}; // no effect
        }

    const auto bound = this->_xc - beta / g;
    const auto u = g > 0. ? bound : this->_xc + this->_r;
    const auto l = g > 0. ? this->_xc - this->_r : bound;

    this->_r = algo::half_nonnegative(u - l);
    this->_xc = l + this->_r;
    return {CUTStatus::success, tsq};
}
