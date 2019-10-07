#include "benchmark/benchmark.h"
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/lmi_oracle.hpp>
#include <ellcpp/oracles/lmi_old_oracle.hpp>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
// #include <xtensor/xarray.hpp>

// using namespace fun;
template <typename Oracle>
class my_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

  private:
    Oracle lmi1;
    Oracle lmi2;
    Arr c;

  public:
    my_oracle(const std::vector<Arr>& F1, Arr B1, const std::vector<Arr>& F2, Arr B2,
        Arr& c)
        : lmi1 {F1, B1}
        , lmi2 {F2, B2}
        , c {c}
    {
    }

    std::tuple<Arr, double, bool> operator()(Arr& x, double t)
    {
        using xt::linalg::dot;

        auto f0 = dot(this->c, x)();
        auto fj1 = f0 - t;
        if (fj1 > 0)
        {
            return std::tuple {this->c, fj1, t};
        }
        auto [g2, fj2] = this->lmi1(x);
        if (g2.shape()[0] > 1 || g2(0) != 0.)
        {
            return std::tuple {std::move(g2), fj2, t};
        }
        auto [g3, fj3] = this->lmi2(x);
        if (g3.shape()[0] > 1 || g3(0) != 0.)
        {
            return std::tuple {std::move(g3), fj3, t};
        }
        return std::tuple {this->c, 0., f0};
    }
};

static void BM_LMI_Lazy(benchmark::State& state)
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

    while (state.KeepRunning())
    {
        auto P = my_oracle<lmi_oracle>(F1, B1, F2, B2, c);
        auto E = ell(10., Arr {0., 0., 0.});
        auto ell_info = cutting_plane_dc(P, E, std::numeric_limits<double>::max());
    }
}

// Register the function as a benchmark
BENCHMARK(BM_LMI_Lazy);

//~~~~~~~~~~~~~~~~

// Define another benchmark
static void BM_LMI_old(benchmark::State& state)
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

    while (state.KeepRunning())
    {
        auto P = my_oracle<lmi_old_oracle>(F1, B1, F2, B2, c);
        auto E = ell(10., Arr {0., 0., 0.});
        auto ell_info = cutting_plane_dc(P, E, std::numeric_limits<double>::max());
    }
}
BENCHMARK(BM_LMI_old);

static void BM_LMI_No_Trick(benchmark::State& state)
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

    while (state.KeepRunning())
    {
        auto P = my_oracle<lmi_oracle>(F1, B1, F2, B2, c);
        auto E = ell(10., Arr {0., 0., 0.});
        E._no_defer_trick = true;
        auto ell_info = cutting_plane_dc(P, E, std::numeric_limits<double>::max());
    }
}

// Register the function as a benchmark
BENCHMARK(BM_LMI_No_Trick);

BENCHMARK_MAIN();
