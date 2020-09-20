#include "benchmark/benchmark.h"
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/lmi_old_oracle.hpp>
#include <ellcpp/oracles/lmi_oracle.hpp>
#include <gsl/span>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>


/*!
 * @brief
 *
 * @tparam Oracle
 */
template <typename Oracle>
class my_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using Cut = std::tuple<Arr, double>;

  private:
    Oracle lmi1;
    Oracle lmi2;
    const Arr c;

  public:
    /*!
     * @brief Construct a new my oracle object
     *
     * @param[in] F1
     * @param[in] B1
     * @param[in] F2
     * @param[in] B2
     * @param[in] c
     */
    my_oracle(gsl::span<const Arr> F1, const Arr& B1, gsl::span<const Arr> F2,
        const Arr& B2, Arr c)
        : lmi1 {F1, B1}
        , lmi2 {F2, B2}
        , c {std::move(c)}
    {
    }

    /*!
     * @brief
     *
     * @param[in] x
     * @param[in] t
     * @return std::tuple<Cut, double>
     */
    std::tuple<Cut, bool> operator()(const Arr& x, double& t)
    {
        const auto f0 = xt::linalg::dot(this->c, x)();
        // const auto f1 = f0 - t;
        if (const auto f1 = f0 - t; f1 > 0)
        {
            return {{this->c, f1}, false};
        }
        const auto cut2 = this->lmi1(x);
        if (cut2)
        {
            return {*cut2, false};
        }
        const auto cut3 = this->lmi2(x);
        if (cut3)
        {
            return {*cut3, false};
        }
        t = f0;
        return {{this->c, 0.}, true};
    }
};

/*!
 * @brief
 *
 * @param[in,out] state
 */
static void BM_LMI_Lazy(benchmark::State& state)
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

    // auto c = Arr {1., -1., 1.};
    const auto F1 = std::vector<Arr> {{{-7., -11.}, {-11., 3.}},
        {{7., -18.}, {-18., 8.}}, {{-2., -8.}, {-8., 1.}}};
    const auto B1 = Arr {{33., -9.}, {-9., 26.}};
    const auto F2 =
        std::vector<Arr> {{{-21., -11., 0.}, {-11., 10., 8.}, {0., 8., 5.}},
            {{0., 10., 16.}, {10., -10., -10.}, {16., -10., 3.}},
            {{-5., 2., -17.}, {2., -6., 8.}, {-17., 8., 6.}}};
    const auto B2 = Arr {{14., 9., 40.}, {9., 91., 10.}, {40., 10., 15.}};

    while (state.KeepRunning())
    {
        auto P = my_oracle<lmi_oracle>(F1, B1, F2, B2, Arr {1., -1., 1.});
        auto E = ell(10., Arr {0., 0., 0.});
        auto t = 1.e100; // std::numeric_limits<double>::max()
        [[maybe_unused]] const auto rslt = cutting_plane_dc(P, E, t);
    }
}

// Register the function as a benchmark
BENCHMARK(BM_LMI_Lazy);

//~~~~~~~~~~~~~~~~

/*!
 * @brief Define another benchmark
 *
 * @param[in,out] state
 */
static void BM_LMI_old(benchmark::State& state)
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

    // auto c = Arr {1., -1., 1.};
    const auto F1 = std::vector<Arr> {{{-7., -11.}, {-11., 3.}},
        {{7., -18.}, {-18., 8.}}, {{-2., -8.}, {-8., 1.}}};
    const auto B1 = Arr {{33., -9.}, {-9., 26.}};
    const auto F2 =
        std::vector<Arr> {{{-21., -11., 0.}, {-11., 10., 8.}, {0., 8., 5.}},
            {{0., 10., 16.}, {10., -10., -10.}, {16., -10., 3.}},
            {{-5., 2., -17.}, {2., -6., 8.}, {-17., 8., 6.}}};
    const auto B2 = Arr {{14., 9., 40.}, {9., 91., 10.}, {40., 10., 15.}};

    while (state.KeepRunning())
    {
        auto P = my_oracle<lmi_old_oracle>(F1, B1, F2, B2, Arr {1., -1., 1.});
        auto E = ell(10., Arr {0., 0., 0.});
        auto t = 1.e100; // std::numeric_limits<double>::max()
        [[maybe_unused]] const auto rslt = cutting_plane_dc(P, E, t);
    }
}
BENCHMARK(BM_LMI_old);

/*!
 * @brief
 *
 * @param[in,out] state
 */
static void BM_LMI_No_Trick(benchmark::State& state)
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

    // const auto c = Arr {1., -1., 1.};
    const auto F1 = std::vector<Arr> {{{-7., -11.}, {-11., 3.}},
        {{7., -18.}, {-18., 8.}}, {{-2., -8.}, {-8., 1.}}};
    const auto B1 = Arr {{33., -9.}, {-9., 26.}};
    const auto F2 =
        std::vector<Arr> {{{-21., -11., 0.}, {-11., 10., 8.}, {0., 8., 5.}},
            {{0., 10., 16.}, {10., -10., -10.}, {16., -10., 3.}},
            {{-5., 2., -17.}, {2., -6., 8.}, {-17., 8., 6.}}};
    const auto B2 = Arr {{14., 9., 40.}, {9., 91., 10.}, {40., 10., 15.}};

    while (state.KeepRunning())
    {
        auto P = my_oracle<lmi_oracle>(F1, B1, F2, B2, Arr {1., -1., 1.});
        auto E = ell(10., Arr {0., 0., 0.});
        E.no_defer_trick = true;
        auto t = 1.e100; // std::numeric_limits<double>::max()
        [[maybe_unused]] const auto rslt = cutting_plane_dc(P, E, t);
    }
}

// Register the function as a benchmark
BENCHMARK(BM_LMI_No_Trick);

BENCHMARK_MAIN();

/*
----------------------------------------------------------
Benchmark                Time             CPU   Iterations
----------------------------------------------------------
BM_LMI_Lazy         131235 ns       131245 ns         4447
BM_LMI_old          196694 ns       196708 ns         3548
BM_LMI_No_Trick     129743 ns       129750 ns         5357
*/