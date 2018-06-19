// -*- coding: utf-8 -*-
#include <catch.hpp>
#include <iostream>
#include <string>
#include <any>

#include <py2cpp/py2cpp.hpp>


TEST_CASE("PY2CPP", "[py2cpp]") {
    // using namespace std::string_literals;

    auto S = py::set{"hello", "world", "test"};
    CHECK(S.contains("test"));
    CHECK(!S.contains("test2"));
    CHECK(py::len(S) == 3);
    for (auto e : S) {
        CHECK(S.contains(e));
    }
    auto& S2 = S;
    auto S3 = S.copy();

    auto M = py::dict<int, std::any>{{1, "hello"}, {8, 5.6}, {9, 4.2}};
    CHECK(M.contains(8));
    CHECK(!M.contains(10));    
    CHECK(py::len(M) == 3);    
    for (auto [e, v] : M) {
        CHECK(M.contains(e));
    }
    CHECK(std::any_cast<double>(M[8]) == 5.6);
    CHECK(std::any_cast<double>(M[9]) != 4.1);
    CHECK(std::any_cast<double>(M.get(2, std::any(-1000.))) == -1000.);

    auto& M2 = M;
    auto M3 = M.copy();

}