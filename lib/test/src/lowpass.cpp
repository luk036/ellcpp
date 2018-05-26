// -*- coding: utf-8 -*-
#include <catch.hpp>
#include <iostream>
#include <tuple>

#include <cmath>
#include <complex>
#include <cutting_plane.hpp>
#include <ell.hpp>
#include <limits>
#include <lowpass_oracle.hpp>

using Arr = xt::xarray<double>;
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
static auto N = 32;
static auto wpass = 0.12 * PI; // end of passband
static auto wstop = 0.20 * PI; // start of stopband
static auto delta0_wpass = 0.125;
static auto delta0_wstop = 0.125;
// maximum passband ripple in dB (+/- around 0 dB)
static auto delta = 20 * std::log10(1 + delta0_wpass);
// stopband attenuation desired in dB
static auto delta2 = 20 * std::log10(delta0_wstop);

// *********************************************************************
// optimization parameters
// *********************************************************************
// rule-of-thumb discretization (from Cheney's Approximation Theory)
static auto m = 15 * N;
static Arr w = xt::linspace<double>(0, PI, m); // omega

// A is the matrix used to compute the power spectrum
// A(w,:) = [1 2*cos(w) 2*cos(2*w) ... 2*cos(N*w)]
static auto An = 2 * xt::cos(xt::linalg::outer(w, xt::arange(1, N)));
static auto A = xt::concatenate(xt::xtuple(xt::ones<double>({m}), An), 1);

// passband 0 <= w <= w_pass
static auto ind_p = xt::where(w <= wpass)[0]; // passband
static double Lp = std::pow(10, -delta / 20);
static double Up = std::pow(10, +delta / 20);
static Arr Ap = xt::view(A, xt::range(0, ind_p.size()), xt::all());

// stopband (w_stop <= w)
static auto ind_s = xt::where(wstop <= w)[0]; // stopband
static double Sp = std::pow(10, delta2 / 20);

using xt::placeholders::_;
static Arr As = xt::view(A, xt::range(ind_s[0], _), xt::all());

// remove redundant contraints
// ind_nr = setdiff(1:m,ind_p)   // fullband less passband
// ind_nr = setdiff(ind_nr, ind_s) // luk: for making parallel cut
// auto ind_nr = np.setdiff1d(xt::arange(m), ind_p);
// auto ind_nr = np.setdiff1d(ind_nr, ind_s);
static auto ind_beg = ind_p[ind_p.size() - 1];
static auto ind_end = ind_s[0];
static Arr Anr = xt::view(A, xt::range(ind_beg, ind_end), xt::all());

static double Lpsq = Lp * Lp;
static double Upsq = Up * Up;
static double Spsq = Sp * Sp;
// ********************************************************************
// optimization
// ********************************************************************

auto run_lowpass(bool use_parallel) {
  Arr r0 = xt::zeros<double>({N}); // initial x0
  // r0[0] = 0;
  auto E = ell(4., r0);
  E._use_parallel = use_parallel;
  auto P = lowpass_oracle(Ap, As, Anr, Lpsq, Upsq);
  auto options = Options();
  options.max_it = 20000;
  options.tol = 1e-4;
  auto [r, Spsq_new, num_iters, feasible, status] =
      cutting_plane_dc(P, E, Spsq, options);
  return std::tuple<bool, unsigned int>{feasible, num_iters};
}

// auto test_lowpass0(benchmark) {
//     result = benchmark(run_lowpass, False)
//     assert result == 13325
// }

// auto test_lowpass1(benchmark) {
//     result = benchmark(run_lowpass, True)
//     assert result == 568
// }

TEST_CASE("Lowpass Filter (w/ parallel cut)", "[lowpass]") {
  auto [feasible, num_iters] = run_lowpass(true);
  CHECK(feasible);
  CHECK(num_iters <= 414);
}

TEST_CASE("Lowpass Filter (w/o parallel cut)", "[lowpass]") {
  auto [feasible, num_iters] = run_lowpass(false);
  CHECK(feasible);
  CHECK(num_iters >= 15110);
}