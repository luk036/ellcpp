#include <catch.hpp>
#include <ellcpp/oracles/chol_ext.hpp>
// #include <xtensor/xarray.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

TEST_CASE("Cholesky test 1", "[chol_ext]")
{
    auto m1 = Arr{{25., 15., -5.}, {15., 18., 0.}, {-5., 0., 11.}};
    auto Q1 = chol_ext(m1.shape()[0]);
    Q1.factorize(m1);
    CHECK(Q1.is_spd());
}

TEST_CASE("Cholesky test 2", "[chol_ext]")
{
    auto m2 = Arr{{18., 22., 54., 42.},
                  {22., -70., 86., 62.},
                  {54., 86., -174., 134.},
                  {42., 62., 134., -106.}};
    auto Q2 = chol_ext(m2.shape()[0]);
    Q2.factorize(m2);
    CHECK(!Q2.is_spd());
    CHECK(Q2.stop == 2);
    // auto ep = Q2.witness();
    // CHECK(ep == 1.);
}

TEST_CASE("Cholesky test 3", "[chol_ext]")
{
    auto m3 = Arr{{0., 15., -5.}, {15., 18., 0.}, {-5., 0., 11.}};
    auto Q3 = chol_ext(m3.shape()[0]);
    Q3.factorize(m3);
    CHECK(!Q3.is_spd());
    auto ep3 = Q3.witness();
    CHECK(Q3.stop == 1);
    // CHECK(v3(0) == 1.);
    CHECK(ep3 == 0.);
}
