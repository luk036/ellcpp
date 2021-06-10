// -*- coding: utf-8 -*-
#include <doctest/doctest.h>
#include <string>

extern auto to_csd(double num, int places = 0) noexcept -> std::string;
extern auto to_csdfixed(double num, unsigned int nnz = 4) noexcept
    -> std::string;
extern auto to_decimal(std::string_view csd_str) noexcept -> double;

using namespace std::string_literals;


TEST_CASE("csd1")
{
    auto csdstr = "+00-00+"s;
    auto csdnumber = to_decimal(csdstr);
    auto csdnew = to_csd(csdnumber);
    // CHECK(csdnew == csdstr);
}


TEST_CASE("csd2")
{
    auto csdstr = "+00-.000+"s;
    auto csdnumber = to_decimal(csdstr);
    auto csdnew = to_csd(csdnumber, 4);
    // CHECK(csdnew == csdstr);
}


TEST_CASE("csd3")
{
    auto csdstr = "+00-.000+"s;
    auto csdnumber = to_decimal(csdstr);
    auto csdnew = to_csdfixed(csdnumber, 3);
    // CHECK(csdnew == csdstr);
}


TEST_CASE("csd4")
{
    auto n = 545.;
    auto csdstr = to_csd(n);
    auto n2 = to_decimal(csdstr);
    CHECK(n == n2);
}
