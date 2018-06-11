/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#define CATCH_CONFIG_MAIN
#include <catch.hpp>
#include <cutting_plane.hpp>
#include <ell.hpp>
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
    int niter, status;
    bool feasible;

    {
        ell E(100.0, Vec{0.0, 0.0});
        profit_oracle P(p, A, alpha, beta, v1, v2, k);
        std::tie(std::ignore, fb, niter, feasible, status) =
            cutting_plane_dc(P, E, 0.0);
        // fmt::print("{:f} {} {} {} \n", fb, niter, feasible, status);
        std::cout << fb << ", " << niter << ", " << feasible << ", " << status
                  << "\n";
        CHECK(niter == 37);
    }

    {
        double ui = 1.0, e1 = 0.003, e2 = 0.007, e3 = 1.0;

        ell E1(100.0, Vec{0.0, 0.0});
        profit_rb_oracle P1(p, A, alpha, beta, v1, v2, k, ui, e1, e2, e3);
        std::tie(std::ignore, fb, niter, feasible, status) =
            cutting_plane_dc(P1, E1, 0.0);
        // fmt::print("{:f} {} {} {} \n", fb, niter, feasible, status);
        std::cout << fb << ", " << niter << ", " << feasible << ", " << status
                  << "\n";
        CHECK(niter == 42);
    }

    {
        ell E2(100.0, Vec{2, 0.0});
        profit_q_oracle P2(p, A, alpha, beta, v1, v2, k);
        std::tie(std::ignore, fb, niter, feasible, status) =
            cutting_plane_q(P2, E2, 0.0);
        // fmt::print("{:f} {} {} {} \n", fb, niter, feasible, status);
        std::cout << fb << ", " << niter << ", " << feasible << ", " << status
                  << "\n";
        CHECK(niter == 28);
    }
}

// TEST_CASE( "Projective Point", "[proj_plane]" ) {
//     CHECK( l.incident({l, m}) );
// }

// int main(int argc, char* argv[]) {
//   //using namespace ModernCppCI;

//   auto result = Catch::Session().run(argc, argv);

//   return (result < 0xff ? result : 0xff);
// }
