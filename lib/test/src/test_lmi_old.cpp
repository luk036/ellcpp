/*
 *  Distributed under the MIT License (See accompanying file /LICENSE )
 */
#include <catch2/catch.hpp>
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/lmi_old_oracle.hpp>
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
    lmi_old_oracle lmi1;
    lmi_old_oracle lmi2;
    Arr c;

  public:
    /**
     * @brief Construct a new my oracle object
     * 
     * @param F1 
     * @param B1 
     * @param F2 
     * @param B2 
     * @param c 
     */
    my_oracle(const M_t& F1, const Arr& B1, const M_t& F2,
        const Arr& B2, const Arr& c)
        : lmi1 {F1, B1}
        , lmi2 {F2, B2}
        , c {c}
    {
    }

    /**
     * @brief 
     * 
     * @param x 
     * @param t 
     * @return std::tuple<Cut, double> 
     */
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

TEST_CASE("LMI (old) test", "[lmi_old_oracle]")
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using M_t = std::vector<Arr>;

    const auto c = Arr {1., -1., 1.};
    const auto F1 = M_t {{{-7., -11.}, {-11., 3.}},
        {{7., -18.}, {-18., 8.}}, {{-2., -8.}, {-8., 1.}}};
    const auto B1 = Arr {{33., -9.}, {-9., 26.}};
    const auto F2 =
        M_t {{{-21., -11., 0.}, {-11., 10., 8.}, {0., 8., 5.}},
            {{0., 10., 16.}, {10., -10., -10.}, {16., -10., 3.}},
            {{-5., 2., -17.}, {2., -6., 8.}, {-17., 8., 6.}}};
    const auto B2 = Arr {{14., 9., 40.}, {9., 91., 10.}, {40., 10., 15.}};

    auto P = my_oracle(F1, B1, F2, B2, c);
    auto E = ell(10., Arr {0., 0., 0.});
    auto [_, ell_info] =
        cutting_plane_dc(P, E, std::numeric_limits<double>::max());
    CHECK(ell_info.feasible);
    CHECK(ell_info.num_iters == 115);
}
