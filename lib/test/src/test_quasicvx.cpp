// -*- coding: utf-8 -*-
#include <cmath>
#include <doctest.h>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <limits>
#include <tuple>

using Arr = xt::xarray<double, xt::layout_type::row_major>;
using Cut = std::tuple<Arr, double>;

/*!
 * @brief
 *
 * @param[in] z
 * @param[in] t
 * @return std::tuple<Cut, double>
 */
std::tuple<Cut, bool> my_quasicvx_oracle(const Arr& z, double& t)
{
    auto sqrtx = z[0];
    auto ly = z[1];

    // constraint 1: exp(x) <= y, or sqrtx**2 <= ly
    auto fj =  sqrtx * sqrtx - ly;
    if (fj > 0.)
    {
        return {{Arr {2*sqrtx, -1.}, fj}, false};
    }

    // constraint 2: x - y >= 1
    auto tmp2 = std::exp(ly);
    auto tmp3 = t * tmp2;
    fj = -sqrtx + tmp3;
    if (fj < 0.)  // feasible
    {
        t = sqrtx / tmp2;
        return {{Arr {-1., sqrtx}, 0}, true};
    }

    return {{Arr {-1., tmp3}, fj}, false};
}

TEST_CASE("Quasiconvex 1, test feasible")
{
    auto E = ell {10., Arr {0., 0.}};
    const auto P = my_quasicvx_oracle;
    auto t = 0.;
    const auto [x, ell_info] = cutting_plane_dc(P, E, t);
    CHECK(ell_info.feasible);
    CHECK(-t == doctest::Approx(-0.4288673397));
    CHECK(x[0]*x[0] == doctest::Approx(0.5029823096));
    CHECK(std::exp(x[1]) == doctest::Approx(1.6536872635));
}

TEST_CASE("Quasiconvex 1, test infeasible 1")
{
    auto E = ell {10., Arr {100., 100.}}; // wrong initial guess
                                          // or ellipsoid is too small
    const auto P = my_quasicvx_oracle;
    auto t = 0.;
    const auto [x, ell_info] = cutting_plane_dc(P, E, t);
    CHECK(!ell_info.feasible);
    CHECK(ell_info.status == CUTStatus::nosoln); // no sol'n
}

TEST_CASE("Quasiconvex 1, test infeasible 2")
{
    auto E = ell {10., Arr {0., 0.}};
    const auto P = my_quasicvx_oracle;
    const auto [x, ell_info] =
        cutting_plane_dc(P, E, 100.); // wrong initial guess
    CHECK(!ell_info.feasible);
}
