#include "profitmaxprob.hpp"
#include <iostream>
#include <yalaa.hpp> // aaf

using ivt = yalaa::details::double_iv_t;
using aaf = yalaa::aff_e_d;

static const double ui = 0.5;
static const double e1 = 0.003, e2 = 0.007;
static const double e3 = 0.5, e4 = 0.5, e5 = 0.5, e6 = 0.5;

/**
 * Profit Maximization Problem, taken from:
 *
 * [1] Aliabadi, Hossein, and Maziar Salahi. "Robust Geometric Programming
 * Approach to Profit Maximization with Interval Uncertainty." Computer
 * Science Journal of Moldova 21.1 (2013): 61.
 */
template <>
profit_max<aaf>::profit_max()
    : gp_base<aaf>()
    , p {ivt(20 - ui * e3, 20 + ui * e3)}
    , alpha {ivt(0.1 - ui * e1, 0.1 + ui * e1)}
    , beta {ivt(0.4 - ui * e2, 0.4 + ui * e2)}
    , v1 {ivt(10.0 - ui * e4, 10.0 + ui * e4)}
    , v2 {ivt(35.0 - ui * e5, 35.0 + ui * e5)}
    , k {ivt(40 - ui * e6, 40.0 + ui * e5)}
    , A {40.0}
    , cd(3)
{
    gp_setup();
}

template <>
profit_max<aaf>::~profit_max()
{
}

int main()
{
    profit_max<aaf> P;
    std::valarray<double> z = {0.0, 1.0, 1.0};
    // z[0] = log(P.obj(z));
    double bf; // for output
    ellipsoid E(z, 200.0);

    STATUS status = ellipsoid_dc(E, P, z, bf, 1000, 1e-6);
    if (status == FOUND)
    {
        std::cout << P.obj(z) << std::endl;
    }
}
