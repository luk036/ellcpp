#include <cmath>
#include <ellcpp/ell1d.hpp>
// #include <tuple>

/* linux-2.6.38.8/include/linux/compiler.h */
#include <cstdio>
#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

/*!
 * @brief
 *
 * @param g
 * @param beta
 * @return ell1d::return_t
 */
ell1d::return_t ell1d::update(const std::tuple<double, double>& cut)
{
    const auto& [g, beta] = cut;

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
