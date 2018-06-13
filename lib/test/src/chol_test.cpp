#include <catch.hpp>
#include <ellcpp/oracles/chol_ext.hpp>
#include <iostream>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

TEST_CASE("Cholesky test", "[chol_ext]") {
    using Arr = xt::xarray<double>;
    using xt::linalg::dot;

    auto m1 = Arr{{25., 15., -5.}, {15., 18., 0.}, {-5., 0., 11.}};
    auto Q1 = chol_ext(m1.shape()[0]);
    Q1.factorize(m1);
    REQUIRE(Q1.is_spd());

    auto m2 = Arr{{18., 22., 54., 42.},
                  {22., -70., 86., 62.},
                  {54., 86., -174., 134.},
                  {42., 62., 134., -106.}};
    auto Q2 = chol_ext(m2.shape()[0]);
    Q2.factorize(m2);
    REQUIRE(!Q2.is_spd());

    auto v = Q2.witness();
    auto p = v.size();

    REQUIRE(p == 2);
}