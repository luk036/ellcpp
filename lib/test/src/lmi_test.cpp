/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include <catch2/catch.hpp>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/lmi_oracle.hpp>
#include <fmt/format.h>
#include <spdlog/sinks/stdout_sinks.h>
#include <spdlog/spdlog.h>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>

/*!
 * @brief my_oracle
 *
 */
class my_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using M_t = std::vector<Arr>;
    using Cut = std::tuple<Arr, double>;

  private:
    lmi_oracle lmi1;
    lmi_oracle lmi2;
    Arr c;

  public:
    /*!
     * @brief Construct a new my oracle object
     *
     * @param F1
     * @param B1
     * @param F2
     * @param B2
     * @param c
     */
    my_oracle(const std::vector<Arr>& F1, const Arr& B1,
        const std::vector<Arr>& F2, const Arr& B2, const Arr& c)
        : lmi1 {F1, B1}
        , lmi2 {F2, B2}
        , c {c}
    {
    }

    /*!
     * @brief
     *
     * @param x
     * @param t
     * @return std::tuple<Cut, double>
     */
    std::tuple<Cut, double> operator()(const Arr& x, double t)
    {
        if (auto cut = this->lmi1(x))
        {
            return {*cut, t};
        }
        if (auto cut = this->lmi2(x))
        {
            return {*cut, t};
        }
        const auto f0 = xt::linalg::dot(this->c, x)();
        const auto f1 = f0 - t;
        if (f1 > 0)
        {
            return {{this->c, f1}, t};
        }
        return {{this->c, 0.}, f0};
    }
};

TEST_CASE("LMI test", "[lmi_oracle]")
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using M_t = std::vector<Arr>;

    const auto c = Arr {1., -1., 1.};
    auto F1 = M_t {{{-7., -11.}, {-11., 3.}}, {{7., -18.}, {-18., 8.}},
        {{-2., -8.}, {-8., 1.}}};
    auto B1 = Arr {{33., -9.}, {-9., 26.}};
    auto F2 = M_t {{{-21., -11., 0.}, {-11., 10., 8.}, {0., 8., 5.}},
        {{0., 10., 16.}, {10., -10., -10.}, {16., -10., 3.}},
        {{-5., 2., -17.}, {2., -6., 8.}, {-17., 8., 6.}}};
    auto B2 = Arr {{14., 9., 40.}, {9., 91., 10.}, {40., 10., 15.}};

    auto P = my_oracle(F1, B1, F2, B2, c);
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

    // create color multi threaded logger
    auto console = spdlog::stdout_color_mt("console");
    auto err_logger = spdlog::stderr_color_mt("stderr");
    spdlog::get("console")->info("loggers can be retrieved from a global "
                                 "registry using the spdlog::get(logger_name)");

    CHECK(ell_info.feasible);
    CHECK(ell_info.num_iters == 113);
}

// TEST_CASE( "Projective Point", "[proj_plane]" ) {
//     REQUIRE( l.incident({l, m}) );
// }

// int main(int argc, char* argv[]) {
//   //using namespace ModernCppCI;

//   auto result = Catch::Session().run(argc, argv);

//   return (result < 0xff ? result : 0xff);
// }