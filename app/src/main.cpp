#include "cutting_plane.hpp"
#include "ell.hpp"
#include "profit_oracle.hpp"
#include "lmi_oracle.hpp"
#include "chol_ext.hpp"
//#include <boost/numeric/ublas/symmetric.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor-blas/xlinalg.hpp>
//#include <fmt/format.h>
#include <iostream>

//#include <boost/numeric/ublas/io.hpp>

// Versions: (latest first)
// http://melpon.org/wandbox/permlink/CcvL0BHfVJhHZH4M
// http://melpon.org/wandbox/permlink/ExQoromITFQ7WOxO
// http://melpon.org/wandbox/permlink/YiLtKIriWtkigZs8

int test1() {
  //namespace bnu = boost::numeric::ublas;
  //using Vec = bnu::vector<double>;
  //using Mat = xt::xarray<double>;
  using Vec = xt::xarray<double>;
 
  double p = 20, A = 40, alpha = 0.1, beta = 0.4;
  double v1 = 10, v2 = 35, k = 30.5;
  double fb;
  int iter, flag, status;

  {
    ell E(100.0, Vec{0.0, 0.0});
    profit_oracle P(p, A, alpha, beta, v1, v2, k);
    std::tie(std::ignore, fb, iter, flag, status) =
        cutting_plane_dc(P, E, 0.0, 200, 1e-4);
    // fmt::print("{:f} {} {} {} \n", fb, iter, flag, status);
    std::cout << fb << ", " << iter << ", " << flag << ", " << status << "\n";
  }

  double ui = 1.0, e1 = 0.003, e2 = 0.007, e3 = 1.0;

  {
    ell E1(100.0, Vec{0.0, 0.0});
    profit_rb_oracle P1(p, A, alpha, beta, v1, v2, k, ui, e1, e2, e3);
    std::tie(std::ignore, fb, iter, flag, status) =
        cutting_plane_dc(P1, E1, 0.0, 200, 1e-4);
    // fmt::print("{:f} {} {} {} \n", fb, iter, flag, status);
    std::cout << fb << ", " << iter << ", " << flag << ", " << status
              << "\n";
  }

  {
    ell E2(100.0, Vec{2, 0.0});
    profit_q_oracle P2(p, A, alpha, beta, v1, v2, k);
    std::tie(std::ignore, fb, iter, flag, status) =
        cutting_plane_q(P2, E2, 0.0, 200, 1e-4);
    // fmt::print("{:f} {} {} {} \n", fb, iter, flag, status);
    std::cout << fb << ", " << iter << ", " << flag << ", " << status
             << "\n";
  }

  return 0;
}

int main() {
  using Arr = xt::xarray<double>;
  using xt::placeholders::_;
  using xt::linalg::dot;

  auto m1 = Arr({{25., 15., -5.},
                {15., 18.,  0.},
                {-5.,  0., 11.}});
  auto Q1 = chol_ext(m1);
  if (!Q1.is_sd()) {
    auto v = Q1.witness();
    auto p = v.size();
    auto sub = xt::range(_, p);
    Arr App = xt::view(m1, sub, sub);
    Arr Appv = dot(App, v);
    auto fj = -dot(v, Appv)();
    std::cout << fj << std::endl;
  }


  auto m2 = Arr({{18., 22.,  54.,  42.},
                {22., -70.,  86.,  62.},
                {54., 86., -174., 134.},
                {42., 62., 134., -106.}});
  auto Q2 = chol_ext(m2);
  if (!Q2.is_sd()) {
    auto v = Q2.witness();
    auto p = v.size();
    auto sub = xt::range(_, p);
    Arr App = xt::view(m2, sub, sub);
    Arr Appv = dot(App, v);
    auto fj = -dot(v, Appv)();
    std::cout << fj << std::endl;
  }

  return 0;
}