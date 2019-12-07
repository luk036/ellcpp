// -*- coding: utf-8 -*-
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <iostream>
#include <tuple>

#include <cmath>
// #include <complex>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/lowpass_oracle.hpp>
#include <ellcpp/utility.hpp>
#include <limits>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xview.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;
// using CArr = xt::xarray<std::complex<double>>;
// using namespace std::literals::complex_literals;

static const auto PI = std::acos(-1);

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
static const auto N            = 32;
static const auto wpass        = 0.12 * PI; // end of passband
static const auto wstop        = 0.20 * PI; // start of stopband
static const auto delta0_wpass = 0.125;
static const auto delta0_wstop = 0.125;
// maximum passband ripple in dB (+/- around 0 dB)
static const auto delta = 20 * std::log10(1 + delta0_wpass);
// stopband attenuation desired in dB
static const auto delta2 = 20 * std::log10(delta0_wstop);

// *********************************************************************
// optimization parameters
// *********************************************************************
// rule-of-thumb discretization (from Cheney's Approximation Theory)
static const auto m = 15 * N;
static const Arr w = xt::linspace<double>(0, PI, m); // omega

// A is the matrix used to compute the power spectrum
// A(w,:) = [1 2*cos(w) 2*cos(2*w) ... 2*cos(N*w)]
static Arr An = 2 * xt::cos(xt::linalg::outer(w, xt::arange(1, N)));
static Arr A  = xt::concatenate(xt::xtuple(xt::ones<double>({m, 1}), An), 1);

// passband 0 <= w <= w_pass
static const auto ind_p = xt::where(w <= wpass)[0]; // passband
static const auto Lp    = std::pow(10, -delta / 20);
static const auto Up    = std::pow(10, +delta / 20);
static const auto Ap    = Arr{xt::view(A, xt::range(0, ind_p.size()), xt::all())};

// stopband (w_stop <= w)
static const auto ind_s = xt::where(wstop <= w)[0]; // stopband
static const auto Sp    = std::pow(10, delta2 / 20);

using xt::placeholders::_;
static Arr As = xt::view(A, xt::range(ind_s[0], _), xt::all());

// remove redundant contraints
// ind_nr = setdiff(1:m,ind_p)   // fullband less passband
// ind_nr = setdiff(ind_nr, ind_s) // luk: for making parallel cut
// const auto ind_nr = np.setdiff1d(xt::arange(m), ind_p);
// const auto ind_nr = np.setdiff1d(ind_nr, ind_s);
static const auto ind_beg = ind_p[ind_p.size() - 1];
static const auto ind_end = ind_s[0];
static const Arr  Anr     = xt::view(A, xt::range(ind_beg, ind_end), xt::all());

static const auto Lpsq = Lp * Lp;
static const auto Upsq = Up * Up;
static const auto Spsq = Sp * Sp;
// ********************************************************************
// optimization
// ********************************************************************

auto run_lowpass(bool use_parallel_cut)
{
    auto r0 = zeros({N}); // initial x0
    // r0[0] = 0;
    auto E       = ell(4., r0);
    auto P       = lowpass_oracle(Ap, As, Anr, Lpsq, Upsq);
    const auto options = Options{20000, 1e-8};
    E._use_parallel_cut = use_parallel_cut;
    auto t = Spsq;
    const auto [r, ell_info] = cutting_plane_dc(P, E, t, options);
    std::cout << r << '\n';
    return std::tuple<bool, unsigned int>{ell_info.feasible, ell_info.num_iters};
}

// auto test_lowpass0(benchmark) {
//     result = benchmark(run_lowpass, False)
//     assert result == 13325
// }

// auto test_lowpass1(benchmark) {
//     result = benchmark(run_lowpass, True)
//     assert result == 568
// }

TEST_CASE("Lowpass Filter (w/ parallel cut)", "[lowpass]")
{
    // void test1() {
    const auto [feasible, num_iters] = run_lowpass(true);
    CHECK(feasible);
    CHECK(num_iters <= 510);
}

TEST_CASE("Lowpass Filter (w/o parallel cut)", "[lowpass]")
{
    // void test2() {
    const auto [feasible, num_iters] = run_lowpass(false);
    CHECK(feasible);
    CHECK(num_iters >= 13333);
}

// int main() {
//     test1();
//     test2();
// }