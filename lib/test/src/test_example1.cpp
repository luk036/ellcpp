// -*- coding: utf-8 -*-
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
std::tuple<Cut, double> my_oracle(const Arr& z, double t)
{
    auto x = z[0];
    auto y = z[1];

    // constraint 1: x + y <= 3
    auto fj = x + y - 3.;
    if (fj > 0.)
    {
        return {{Arr {1., 1.}, fj}, t};
    }

    // constraint 2: x - y >= 1
    fj = -x + y + 1.;
    if (fj > 0.)
    {
        return {{Arr {-1., 1.}, fj}, t};
    }

    // objective: maximize x + y
    auto f0 = x + y;
    fj = t - f0;
    if (fj < 0.)
    {
        fj = 0.;
        t = f0;
    }
    return {{Arr {-1., -1.}, fj}, t};
}

TEST_CASE("Example 1, test feasible")
{
    auto E = ell {10., Arr {0., 0.}};
    const auto P = my_oracle;
    auto t = -1.e100; // std::numeric_limits<double>::min()
    const auto [x, ell_info] = cutting_plane_dc(P, E, t);
    CHECK(x[0] >= 0.);
    CHECK(ell_info.feasible);
}

TEST_CASE("Example 1, test infeasible 1")
{
    auto E = ell {10., Arr {100., 100.}}; // wrong initial guess
                                          // or ellipsoid is too small
    const auto P = my_oracle;
    auto t = -1.e100; // std::numeric_limits<double>::min()
    const auto [x, ell_info] = cutting_plane_dc(P, E, t);
    CHECK(!ell_info.feasible);
    CHECK(ell_info.status == CUTStatus::nosoln); // no sol'n
}

TEST_CASE("Example 1, test infeasible 2")
{
    auto E = ell {10., Arr {0., 0.}};
    const auto P = my_oracle;
    const auto [x, ell_info] =
        cutting_plane_dc(P, E, 100); // wrong initial guess
    CHECK(!ell_info.feasible);
}
