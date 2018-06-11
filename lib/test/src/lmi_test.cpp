/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include <catch.hpp>
#include <iostream>
#include <tuple>

#include <cutting_plane.hpp>
#include <ell.hpp>
#include <lmi_oracle.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

// using namespace fun;
class my_oracle {
    using Arr = xt::xarray<double>;

  private:
    Arr _c;
    lmi_oracle _lmi1;
    lmi_oracle _lmi2;

  public:
    explicit my_oracle(Arr &c, Arr &F1, Arr &B1, Arr &F2, Arr &B2)
        : _c{c}, _lmi1{F1, B1}, _lmi2{F2, B2} {}

    auto operator()(Arr &x, double t) {
        using xt::linalg::dot;

        auto f0 = dot(_c, x)();
        auto fj1 = f0 - t;
        if (fj1 > 0.0) {
            return std::tuple{_c, fj1, t};
        }

        auto [g2, fj2] = _lmi1.chk_spd(x);
        if (fj2 > 0.) {
            return std::tuple{g2, fj2, t};
        }

        auto [g3, fj3] = _lmi2.chk_spd(x);
        if (fj3 > 0.) {
            return std::tuple{g3, fj3, t};
        }

        return std::tuple{_c, 0.0, f0};
    }
};

TEST_CASE("LMI test", "[lmi_oracle]") {
    using Arr = xt::xarray<double>;

    auto c = Arr{1., -1., 1.};
    auto F1 = Arr{{{-7., -11.}, {-11., 3.}},
                  {{7., -18.}, {-18., 8.}},
                  {{-2., -8.}, {-8., 1.}}};
    auto B1 = Arr{{33., -9.}, {-9., 26.}};
    auto F2 = Arr{{{-21., -11., 0.}, {-11., 10., 8.}, {0., 8., 5.}},
                  {{0., 10., 16.}, {10., -10., -10.}, {16., -10., 3.}},
                  {{-5., 2., -17.}, {2., -6., 8.}, {-17., 8., 6.}}};
    auto B2 = Arr{{14., 9., 40.}, {9., 91., 10.}, {40., 10., 15.}};

    auto P = my_oracle(c, F1, B1, F2, B2);
    auto E = ell(10.0, Arr{0.0, 0.0, 0.0});

    // double fb;
    // int niter, feasible, status;
    // Arr xb;

    auto [xb, fb, niter, feasible, status] = cutting_plane_dc(P, E, 100.0);
    // fmt::print("{:f} {} {} {} \n", fb, niter, feasible, status);
    std::cout << xb << "\n";
    std::cout << fb << ", " << niter << ", " << feasible << ", " << status
              << "\n";

    REQUIRE(feasible);
    REQUIRE(niter == 115);
}

// TEST_CASE( "Projective Point", "[proj_plane]" ) {
//     REQUIRE( l.incident({l, m}) );
// }

// int main(int argc, char* argv[]) {
//   //using namespace ModernCppCI;

//   auto result = Catch::Session().run(argc, argv);

//   return (result < 0xff ? result : 0xff);
// }