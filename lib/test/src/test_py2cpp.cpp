// -*- coding: utf-8 -*-
#include <catch.hpp>
#include <iostream>
#include <string>

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

    auto M = py::dict<int, double>{{1, 4.5}, {8, 5.6}, {9, 4.2}};
    CHECK(M.contains(8));
    CHECK(!M.contains(10));    
    CHECK(py::len(M) == 3);    
    for (auto [e, v] : M) {
        CHECK(M.contains(e));
    }
    CHECK(M[8] == 5.6);
    CHECK(M.get(2, -1000.) == -1000.);
}