/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include <catch2/catch.hpp>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/lmi_oracle.hpp>
#include <fmt/format.h>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
// #include <xtensor/xarray.hpp>

// using namespace fun;
class my_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Cut = std::tuple<Arr, double>;

  private:
    lmi_oracle lmi1;
    lmi_oracle lmi2;
    Arr c;

  public:
    // my_oracle(const std::vector<Arr> &F1, Arr &B1,
    //           const std::vector<Arr> &F2, Arr &B2, Arr &c)
    //     : lmi1{F1, B1}, lmi2{F2, B2}, c{c} {}

    my_oracle(std::vector<Arr>&& F1, Arr&& B1, std::vector<Arr>&& F2, Arr&& B2,
        Arr& c)
        : lmi1 {std::forward<std::vector<Arr>>(F1), std::forward<Arr>(B1)}
        , lmi2 {std::forward<std::vector<Arr>>(F2), std::forward<Arr>(B2)}
        , c {c}
    {
    }

    std::tuple<Cut, double> operator()(const Arr& x, double t)
    {
        using xt::linalg::dot;

        auto f0 = dot(this->c, x)();
        auto fj1 = f0 - t;
        if (fj1 > 0)
        {
            return {{this->c, fj1}, t};
        }
        if (auto cut = this->lmi1(x))
        {
            return {*cut, t};
        }
        if (auto cut = this->lmi2(x))
        {
            return {*cut, t};
        }
        return {{this->c, 0.}, f0};
    }
};

TEST_CASE("LMI test", "[lmi_oracle]")
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

    auto c = Arr {1., -1., 1.};
    auto F1 = std::vector<Arr> {{{-7., -11.}, {-11., 3.}},
        {{7., -18.}, {-18., 8.}}, {{-2., -8.}, {-8., 1.}}};
    auto B1 = Arr {{33., -9.}, {-9., 26.}};
    auto F2 =
        std::vector<Arr> {{{-21., -11., 0.}, {-11., 10., 8.}, {0., 8., 5.}},
            {{0., 10., 16.}, {10., -10., -10.}, {16., -10., 3.}},
            {{-5., 2., -17.}, {2., -6., 8.}, {-17., 8., 6.}}};
    auto B2 = Arr {{14., 9., 40.}, {9., 91., 10.}, {40., 10., 15.}};

    auto P = my_oracle(
        std::move(F1), std::move(B1), std::move(F2), std::move(B2), c);
    auto E = ell(10., Arr {0., 0., 0.});

    // double fb;
    // int niter, feasible, status;
    // Arr xb;

    auto [_, ell_info] =
        cutting_plane_dc(P, E, std::numeric_limits<double>::max());
    fmt::print(
        "{:f} {} {} \n", ell_info.value, ell_info.num_iters, ell_info.feasible);
    // std::cout << "LMI xbest: " << xb << "\n";
    // std::cout << "LMI result: " << fb << ", " << niter << ", " << feasible <<
    // ", " << status
    //           << "\n";

    CHECK(ell_info.feasible);
    CHECK(ell_info.num_iters == 115);
}

// TEST_CASE( "Projective Point", "[proj_plane]" ) {
//     REQUIRE( l.incident({l, m}) );
// }

// int main(int argc, char* argv[]) {
//   //using namespace ModernCppCI;

//   auto result = Catch::Session().run(argc, argv);

//   return (result < 0xff ? result : 0xff);
// }