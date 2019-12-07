#include <cmath>
#include <ellcpp/ell1d.hpp>
// #include <tuple>

/*!
 * @brief
 *
 * @param cut
 * @return ell1d::return_t
 */
ell1d::return_t ell1d::update(const std::tuple<double, double>& cut)
{
    const auto& [g, beta] = cut;

    const auto tau = std::abs(this->_r * g);
    const auto tsq = tau * tau;

    [[unlikely]] if (beta == 0.)
    {
        this->_r /= 2;
        this->_xc += g > 0. ? -this->_r : this->_r;
        return {0, tsq};
    }
    [[unlikely]] if (beta > tau)
    {
        return {1, tsq}; // no sol'n
    }
    [[unlikely]] if (beta < -tau)
    {
        return {3, tsq}; // no effect
    }

    const auto bound = this->_xc - beta / g;
    const auto u = g > 0. ? bound : this->_xc + this->_r;
    const auto l = g > 0. ? this->_xc - this->_r : bound;

    this->_r = (u - l) / 2;
    this->_xc = l + this->_r;
    return {0, tsq};
}
