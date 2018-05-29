// -*- coding: utf-8 -*-
#include <iostream>
#include <tuple>
#include <cutting_plane.hpp>
#include <ell.hpp>
#include <qmi_oracle.hpp>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xmath.hpp>
// #include <vector>

using Arr = xt::xarray<double>;

auto construct_distance_matrix(const Arr &s) {
    auto n = s.shape()[0];
    // c = cvx.Variable(m)
    Arr D = xt::zeros<double>({n, n});
    for (auto i=0u; i<n; ++i) {
        for (auto j=i+1; i<n; ++i) {
            Arr h = xt::view(s, j, xt::all()) - xt::view(s, i, xt::all());
            auto d = xt::sqrt(xt::linalg::dot(h, h))();
            D(i, j) = d;
            D(j, i) = d;
        }
    }
    return D;
}

Arr lsq_corr_poly(const Arr &Y, const Arr &s, std::size_t m) {
    auto n = s.shape()[0];
    Arr a = xt::zeros<double>({m});
    Arr D1 = construct_distance_matrix(s);

    Arr D = xt::ones<double>({n, n});
    Arr Sig = xt::zeros<double>({m, n, n});
    //Sig.push_back(D);
    xt::view(Sig, 0, xt::all(), xt::all()) = D;

    for (auto i=0u; i<m-1; ++i) {
        D *= D1;
        // Sig.push_back(D);
        xt::view(Sig, i, xt::all(), xt::all()) = D;
    }
    // Sig.reverse();

    // P = mtx_norm_oracle(Sig, Y, a)
    auto Q = qmi_oracle(Sig, Y);
    auto E = ell(10., a);
    auto P = bsearch_adaptor(Q, E);
    // double normY = xt::norm_l2(Y);
    auto [fb, niter, feasible] = bsearch(P, std::tuple{0., 100.*100.});

    std::cout << niter << ", " << feasible << '\n';
    a = P.x_best();
    return a;
//  return prob.is_dcp()
}
