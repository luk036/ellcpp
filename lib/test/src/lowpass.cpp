// -*- coding: utf-8 -*-
#include <catch.hpp>
#include <iostream>
#include <tuple>

#include <cmath>
#include <complex>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/lowpass_oracle.hpp>
#include <limits>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xview.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;
// using CArr = xt::xarray<std::complex<double>>;
// using namespace std::literals::complex_literals;

static double PI = std::acos(-1);

// Modified from CVX code by Almir Mutapcic in 2006.
// Adapted in 2010 for impulse response peak-minimization by convex iteration by
// Christine Law.
//
// "FIR Filter Design via Spectral Factorization and Convex Optimization"
// by S.-P. Wu, S. Boyd, and L. Vandenberghe
//
// Designs an FIR lowpass filter using spectral factorization method with
// constraint on maximum passband ripple and stopband attenuation:
//
//   minimize   max |H(w)|                      for w in stopband
//       s.t.   1/delta <= |H(w)| <= delta      for w in passband
//
// We change variables via spectral factorization method and get:
//
//   minimize   max R(w)                          for w in stopband
//       s.t.   (1/delta)**2 <= R(w) <= delta**2  for w in passband
//              R(w) >= 0                         for all w
//
// where R(w) is squared magnitude frequency response
// (and Fourier transform of autocorrelation coefficients r).
// Variables are coeffients r and G = hh' where h is impulse response.
// delta is allowed passband ripple.
// This is a convex problem (can be formulated as an SDP after sampling).

// rand('twister',sum(100*clock))
// randn('state',sum(100*clock))

// *********************************************************************
// filter specs (for a low-pass filter)
// *********************************************************************
// number of FIR coefficients (including zeroth)
static const int N = 48;
static const double wpass = 0.12 * PI; // end of passband
static const double wstop = 0.20 * PI; // start of stopband
static const double delta0_wpass = 0.125;
static const double delta0_wstop = 0.125;
// maximum passband ripple in dB (+/- around 0 dB)
static const double delta = 20 * std::log10(1 + delta0_wpass);
// stopband attenuation desired in dB
static const double delta2 = 20 * std::log10(delta0_wstop);

// *********************************************************************
// optimization parameters
// *********************************************************************
// rule-of-thumb discretization (from Cheney's Approximation Theory)
static const int m = 15 * N;
static Arr w = xt::linspace<double>(0, PI, m); // omega

// A is the matrix used to compute the power spectrum
// A(w,:) = [1 2*cos(w) 2*cos(2*w) ... 2*cos(N*w)]
static Arr An = 2 * xt::cos(xt::linalg::outer(w, xt::arange(1, N)));
static Arr A = xt::concatenate(xt::xtuple(xt::ones<double>({m, 1}), An), 1);

// passband 0 <= w <= w_pass
static auto ind_p = xt::where(w <= wpass)[0]; // passband
static const double Lp = std::pow(10, -delta / 20);
static const double Up = std::pow(10, +delta / 20);
static Arr Ap = xt::view(A, xt::range(0, ind_p.size()), xt::all());

// stopband (w_stop <= w)
static auto ind_s = xt::where(wstop <= w)[0]; // stopband
static const double Sp = std::pow(10, delta2 / 20);

using xt::placeholders::_;
static Arr As = xt::view(A, xt::range(ind_s[0], _), xt::all());

// remove redundant contraints
// ind_nr = setdiff(1:m,ind_p)   // fullband less passband
// ind_nr = setdiff(ind_nr, ind_s) // luk: for making parallel cut
// auto ind_nr = np.setdiff1d(xt::arange(m), ind_p);
// auto ind_nr = np.setdiff1d(ind_nr, ind_s);
static auto ind_beg = ind_p[ind_p.size() - 1];
static auto ind_end = ind_s[0];
static Arr Anr = xt::view(A, xt::range(ind_beg + 1, ind_end), xt::all());

static const double Lpsq = Lp * Lp;
static const double Upsq = Up * Up;
static const double Spsq = Sp * Sp;
// ********************************************************************
// optimization
// ********************************************************************

auto run_lowpass(bool use_parallel_cut) {
    auto r0 = Arr{xt::zeros<double>({N})}; // initial x0
    auto E = ell(40., r0);
    auto P = lowpass_oracle(Ap, As, Anr, Lpsq, Upsq);
    auto options = Options();

    options.max_it = 50000;
    E._use_parallel_cut = use_parallel_cut;
    // options.tol = 1e-8;

    auto [r, Spsq_new, feasible, num_iters, status] =
        cutting_plane_dc(P, E, Spsq, options);
    // std::cout << "lowpass r: " << r << '\n';
    // auto Ustop = 20 * std::log10(std::sqrt(Spsq_new));
    // std::cout << "Min attenuation in the stopband is " << Ustop << " dB.\n";
    return std::tuple{feasible, num_iters};
}

TEST_CASE("Lowpass Filter (w/ parallel cut)", "[lowpass]") {
    auto [feasible, num_iters] = run_lowpass(true);
    CHECK(feasible);
    CHECK(num_iters <= 634);
}

TEST_CASE("Lowpass Filter (w/o parallel cut)", "[lowpass]") {
    auto [feasible, num_iters] = run_lowpass(false);
    CHECK(feasible);
    CHECK(num_iters >= 7479);
}
