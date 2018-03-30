/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include <cutting_plane.hpp>
#include <ell.hpp>
#include <lmi_oracle.hpp>
#include <profit_oracle.hpp>
//#include <boost/numeric/ublas/symmetric.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

// using namespace fun;

TEST_CASE("Profit Test 1", "[profit]") {
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
    REQUIRE(iter == 37);
  }

  double ui = 1.0, e1 = 0.003, e2 = 0.007, e3 = 1.0;

  {
    ell E1(100.0, Vec{0.0, 0.0});
    profit_rb_oracle P1(p, A, alpha, beta, v1, v2, k, ui, e1, e2, e3);
    std::tie(std::ignore, fb, iter, flag, status) =
        cutting_plane_dc(P1, E1, 0.0, 200, 1e-4);
    // fmt::print("{:f} {} {} {} \n", fb, iter, flag, status);
    std::cout << fb << ", " << iter << ", " << flag << ", " << status << "\n";
    REQUIRE(iter == 42);
  }

  {
    ell E2(100.0, Vec{2, 0.0});
    profit_q_oracle P2(p, A, alpha, beta, v1, v2, k);
    std::tie(std::ignore, fb, iter, flag, status) =
        cutting_plane_q(P2, E2, 0.0, 200, 1e-4);
    // fmt::print("{:f} {} {} {} \n", fb, iter, flag, status);
    std::cout << fb << ", " << iter << ", " << flag << ", " << status << "\n";
    REQUIRE(iter == 28);
  }
}

// TEST_CASE( "Projective Point", "[proj_plane]" ) {
//     REQUIRE( l.incident({l, m}) );
// }

// int main(int argc, char* argv[]) {
//   //using namespace ModernCppCI;

//   auto result = Catch::Session().run(argc, argv);

//   return (result < 0xff ? result : 0xff);
// }