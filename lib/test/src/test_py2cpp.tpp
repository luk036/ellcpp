// -*- coding: utf-8 -*-
#include <boost/any.hpp>
#include <doctest/doctest.h>
// #include <iostream>
#include <string>

#include <py2cpp/py2cpp.hpp>

TEST_CASE("PY2CPP")
{
    // using namespace std::string_literals;

    py::set S{"hello", "world", "test"};
    CHECK(S.contains("test"));
    CHECK(!S.contains("test2"));
    CHECK(py::len(S) == 3);
    for (const auto& e : S)
    {
        CHECK(S.contains(e));
    }
    auto S3 = S.copy();

    py::dict<int, boost::any> M{{1, "hello"}, {8, 5.6}, {9, 4.2}};
    CHECK(M.contains(8));
    CHECK(!M.contains(10));
    CHECK(py::len(M) == 3);
    for (auto [e, _] : M.items())
    {
        CHECK(M.contains(e));
    }
    for (const auto& e : M)
    {
        CHECK(M.contains(e));
    }
    CHECK(boost::any_cast<double>(M[8]) == 5.6);
    CHECK(boost::any_cast<double>(M[9]) != 4.1);
    CHECK(boost::any_cast<double>(M.get(2, boost::any(-1000.))) == -1000.);

    auto M3 = M.copy();
}