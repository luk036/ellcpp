// -*- coding: utf-8 -*-
#include <catch.hpp>
#include <iostream>
#include <tuple>

#include <complex>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <limits>

using Arr = xt::xarray<double>;
// using CArr = xt::xarray<std::complex<double>>;
// using namespace std::literals::complex_literals;

static double PI = std::acos(-1);

// ********************************************************************
// Problem specs.
// ********************************************************************
// Number of FIR coefficients (including the zeroth one).
static auto n = 20;

// Rule-of-thumb frequency discretization (Cheney's Approx. Theory book).
static auto m = 15 * n;
static Arr w = xt::linspace<double>(0, PI, m);

// ********************************************************************
// Construct the desired filter.
// ********************************************************************
// Fractional delay.
static auto D = 8.25; // Delay value.
// CArr Hdes = xt::exp(D * w); // Desired frequency response.

// Gaussian filter with linear phase. (Uncomment lines below for this design.)
// var = 0.05
// Hdes = 1/(xt::sqrt(2*PI*var)) * xt::exp(-xt::square(w-PI/2)/(2*var))
// Hdes = xt::multiply(Hdes, xt::exp(-1j*n/2*w))

// A is the matrix used to compute the frequency response
// from a vector of filter coefficients:
//     A[w,:] = [1 exp(-j*w) exp(-j*2*w) ... exp(-j*n*w)]
// CArr A = xt::exp(xt::linalg::kron(w, xt::arange(n)));

// Presently CVXPY does not do complex-valued math, so the
// problem must be formatted into a real-valued representation.

// Split Hdes into a real part, and an imaginary part.
// Arr Hdes_r = xt::real(Hdes);
// Arr Hdes_i = xt::imag(Hdes);
static Arr Hdes_theta = D * w;
static Arr Hdes_r = xt::cos(Hdes_theta);
static Arr Hdes_i = -xt::sin(Hdes_theta);
static Arr A_theta = xt::linalg::outer(w, xt::arange(n));

// Split A into a real part, and an imaginary part.
// Arr A_R = xt::real(A);
// Arr A_I = xt::imag(A);
static Arr A_R = xt::cos(A_theta);
static Arr A_I = -xt::sin(A_theta);

// Optimal Chebyshev filter formulation.
class my_fir_oracle {
  public:
    auto operator()(const Arr &h, double t) const {
        auto fmax = std::numeric_limits<double>::min();
        // auto imax = -1;
        Arr gmax = xt::zeros<double>({n});

        for (auto i = 0u; i < m; ++i) {
            auto a_R = xt::view(A_R, i, xt::all());
            auto a_I = xt::view(A_I, i, xt::all());
            double H_r = Hdes_r[i];
            double H_i = Hdes_i[i];
            double t_r = xt::linalg::dot(a_R, h)() - H_r;
            double t_i = xt::linalg::dot(a_I, h)() - H_i;
            double fj = t_r * t_r + t_i * t_i;
            if (fj >= t) {
                Arr g = 2 * (t_r * a_R + t_i * a_I);
                return std::tuple{std::move(g), fj - t, t};
            }
            if (fmax < fj) {
                fmax = fj;
                // imax = i;
                gmax = 2 * (t_r * a_R + t_i * a_I);
            }
        }

        return std::tuple{std::move(gmax), 0., fmax};
    }
};

// def test_firfilter() {
TEST_CASE("FIR Filter", "[firfilter]") {

    Arr h0 = xt::zeros<double>({n}); // initial x0
    auto E = ell(40., h0);
    auto P = my_fir_oracle();
    auto [hb, fb, niter, feasible, status] = cutting_plane_dc(P, E, 100.);

    CHECK(feasible);
    std::cout << "optimal value " << fb << "\n";
    std::cout << "optimal sol'n " << hb << "\n";

    // fmt = '{ {f} {} {} {}'
    // print(prob1.optim_var)
    // print(fmt.format(prob1.optim_vale, prob1.solver_stats.num_iters))

    // print 'Problem status:', feasible
    // if feasible != 1:
    //    raise Exception('ELL Error')

    // hv = xt::asmatrix(prob1.optim_var).T
}
