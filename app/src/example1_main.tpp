// -*- coding: utf-8 -*-
#include <iostream>
#include <tuple>

#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

auto my_oracle(const Arr &z, double t) {
  auto x = z[0], y = z[1];

  // constraint 1: x + y <= 3
  auto fj = x + y - 3;
  if (fj > 0) {
    return std::tuple{Arr{1., 1.}, fj, t};
  }
  // constraint 2: x - y >= 1
  fj = -x + y + 1;
  if (fj > 0) {
    return std::tuple{Arr{-1., 1.}, fj, t};
  }
  // objective: maximize x + y
  auto f0 = x + y;
  fj = t - f0;
  if (fj < 0) {
    fj = 0.;
    t = f0;
  }
  return std::tuple{Arr{-1., -1.}, fj, t};
}

int main() {
  auto x0 = Arr{0., 0.}; // initial x0
  auto E = ell(10., x0);
  auto P = my_oracle;
  [[maybe_unused]] const auto [_, ell_info] = cutting_plane_dc(P, E, -100.);
  std::cout << ell_info.num_iters << "," << ell_info.feasible << "," << ell_info.status << "\n";
  std::cout << ell_info.val << "\n";
}