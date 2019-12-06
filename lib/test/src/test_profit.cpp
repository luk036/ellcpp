/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include <catch2/catch.hpp>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/profit_oracle.hpp>
#include <xtensor/xarray.hpp>

// using namespace fun;

TEST_CASE("Profit Test 1", "[profit]")
{
    using Vec = xt::xarray<double, xt::layout_type::row_major>;

    const auto p = 20.;
    const auto A = 40.;
    const auto k = 30.5;
    const auto a = Vec {0.1, 0.4};
    const auto v = Vec {10., 35.};

    {
        // auto E = ell(100., Vec {0., 0.});
        auto P = profit_oracle(p, A, k, a, v);
        [[maybe_unused]] const auto [_, ell_info] =
            cutting_plane_dc(P, ell {100., Vec {0., 0.}}, 0.);
        CHECK(ell_info.num_iters == 37);
    }

    {
        auto P = profit_rb_oracle(p, A, k, a, v, Vec {0.003, 0.007}, 1.);
        [[maybe_unused]] const auto [_, ell_info] =
            cutting_plane_dc(P, ell {100., Vec {0., 0.}}, 0.);
        CHECK(ell_info.num_iters == 42);
    }

    {
        auto P = profit_q_oracle(p, A, k, a, v);
        [[maybe_unused]] const auto [_, ell_info] =
            cutting_plane_q(P, ell(100., Vec {2, 0.}), 0.);
        CHECK(ell_info.num_iters == 28);
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
