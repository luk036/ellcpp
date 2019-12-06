// -*- coding: utf-8 -*-
#include <catch2/catch.hpp>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/utility.hpp>
#include <limits>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xview.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

static double PI = std::acos(-1);

// ********************************************************************
// Problem specs.
// ********************************************************************
// Number of FIR coefficients (including the zeroth one).
static auto n = 20U;

// Rule-of-thumb frequency discretization (Cheney's Approx. Theory book).
static auto m = 15U * n;

// ********************************************************************
// Construct the desired filter.
// ********************************************************************
// Fractional delay.
static auto D = 8.25; // Delay value.
// CArr Hdes = xt::exp(D * w); // Desired frequency response.

// Optimal Chebyshev filter formulation.
class my_fir_oracle
{
    using Cut = std::tuple<Arr, double>;

    // Gaussian filter with linear phase. (Uncomment lines below for this
    // design.) var = 0.05 Hdes = 1/(xt::sqrt(2 * PI * var)) *
    // xt::exp(-xt::square(w-PI/2)/(2 * var)) Hdes = xt::multiply(Hdes,
    // xt::exp(-1j * n/2 * w))
    Arr w = xt::linspace<double>(0, PI, m);

    // A is the matrix used to compute the frequency response
    // from a vector of filter coefficients:
    //     A[w,:] = [1 exp(-j * w) exp(-j * 2 * w) ... exp(-j * n * w)]
    // CArr A = xt::exp(xt::linalg::kron(w, xt::arange(n)));

    // Presently CVXPY does not do complex-valued math, so the
    // problem must be formatted into a real-valued representation.

    // Split Hdes into a real part, and an imaginary part.
    // Arr Hdes_r = xt::real(Hdes);
    // Arr Hdes_i = xt::imag(Hdes);
    Arr Hdes_theta = D * w;
    Arr Hdes_r = xt::cos(Hdes_theta);
    Arr Hdes_i = -xt::sin(Hdes_theta);
    Arr A_theta = xt::linalg::outer(w, xt::arange(n));

    // Split A into a real part, and an imaginary part.
    // Arr A_R = xt::real(A);
    // Arr A_I = xt::imag(A);
    Arr A_R = xt::cos(A_theta);
    Arr A_I = -xt::sin(A_theta);

  public:
    auto operator()(const Arr& h, double t) const -> std::tuple<Cut, double>
    {
        auto fmax = std::numeric_limits<double>::min();
        auto gmax = zeros({n});

        for (auto i = 0U; i < m; ++i)
        {
            auto a_R = Arr {xt::view(A_R, i, xt::all())};
            auto a_I = Arr {xt::view(A_I, i, xt::all())};
            auto H_r = Hdes_r[i];
            auto H_i = Hdes_i[i];
            auto t_r = xt::linalg::dot(a_R, h)() - H_r;
            auto t_i = xt::linalg::dot(a_I, h)() - H_i;
            auto fj = t_r * t_r + t_i * t_i;
            if (fj >= t)
            {
                auto g = Arr {2. * (t_r * a_R + t_i * a_I)};
                return {{std::move(g), fj - t}, t};
            }
            if (fmax < fj)
            {
                fmax = fj;
                // imax = i;
                gmax = 2. * (t_r * a_R + t_i * a_I);
            }
        }

        return {{std::move(gmax), 0.}, fmax};
    }
};

TEST_CASE("FIR Filter", "[firfilter]")
{
    auto E = ell(40., zeros({n}));
    auto P = my_fir_oracle();
    auto t = std::numeric_limits<double>::max();
    const auto [_, ell_info] = cutting_plane_dc(P, E, t);
    CHECK(ell_info.feasible);
}
