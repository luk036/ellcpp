// -*- coding: utf-8 -*-
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/qmi_oracle.hpp>
#include <iostream>
#include <tuple>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xnorm.hpp>
#include <vector>

using Arr = xt::xarray<double>;

/**
 * @brief 
 * 
 * @param s 
 * @return Arr 
 */
Arr construct_distance_matrix(const Arr &s) {
    auto n = s.shape()[0];
    // c = cvx.Variable(m)
    Arr D = xt::zeros<double>({n, n});
    for (auto i = 0u; i < n; ++i) {
        for (auto j = i + 1; i < n; ++i) {
            Arr h = xt::view(s, j, xt::all()) - xt::view(s, i, xt::all());
            double d = xt::sqrt(xt::linalg::dot(h, h))();
            D(i, j) = d;
            D(j, i) = d;
        }
    }
    return D;
}

class basis_oracle2 {
    using Arr = xt::xarray<double>;
    using shape_type = Arr::shape_type;

  private:
    qmi_oracle _qmi;

  public:
    explicit basis_oracle2(const std::vector<Arr> &F, const Arr &F0)
        : _qmi(F, F0) {}

    auto operator()(const Arr &x, double t) {
        auto n = x.shape()[0];
        Arr g = xt::zeros<double>({n});
        g(n-1) = 1.;
        auto tc = x(n-1);
        auto f0 = tc - t;
        if (f0 > 0.) {
            return std::tuple{g, f0, t};
        }

        this->_qmi.update(tc);
        auto [g1, fj, feasible] = this->_qmi(xt::view(x, xt::range(0, n-1)));
        if (!feasible) {
            xt::view(g, xt::range(0, n-1)) = g1;
            auto v = this->_qmi._Q.witness();
            g(n-1) = -xt::linalg::dot(v,v)();
            return std::tuple{g, fj, t};
        }

        return std::tuple{g, 0., tc};
    }
};


auto lsq_corr_core2(const Arr &Y, std::size_t m, basis_oracle2 &P) {
    auto normY = 100.;
    auto normY2 = 32*normY*normY;
    Arr val = 256*xt::ones<double>({m + 1});
    val(m) = normY2*normY2;
    Arr x = xt::zeros<double>({m + 1});
    x(m) = normY2/2;
    auto E = ell(val, x);
    auto [x_best, fb, num_iters, feasible, status] = cutting_plane_dc(P, E, normY2);
    assert(feasible);
    Arr a = xt::view(x_best, xt::range(0, m));
    return std::tuple{num_iters, a};
}


std::size_t lsq_corr_poly2(const Arr &Y, const Arr &s, std::size_t m) {
    auto n = s.shape()[0];
    Arr D1 = construct_distance_matrix(s);

    Arr D = xt::ones<double>({n, n});
    // Arr Sig = xt::zeros<double>({m, n, n});
    // xt::view(Sig, 0, xt::all(), xt::all()) = D;
    std::vector<Arr> Sig;
    Sig.reserve(m);
    Sig.push_back(D);

    for (auto i = 0u; i < m - 1; ++i) {
        D *= D1;
        Sig.push_back(D);
        // xt::view(Sig, i, xt::all(), xt::all()) = D;
    }
    // Sig.reverse();

    auto P = basis_oracle2(Sig, Y);
    auto [num_iters, a] = lsq_corr_core2(Y, m, P);
    std::cout << a << "\n";
    return num_iters;
}

std::size_t lsq_corr_poly(const Arr &Y, const Arr &s, std::size_t m) {
    auto n = s.shape()[0];
    Arr D1 = construct_distance_matrix(s);

    Arr D = xt::ones<double>({n, n});
    // Arr Sig = xt::zeros<double>({m, n, n});
    // xt::view(Sig, 0, xt::all(), xt::all()) = D;
    std::vector<Arr> Sig;
    Sig.reserve(m);
    Sig.push_back(D);

    for (auto i = 0u; i < m - 1; ++i) {
        D *= D1;
        // xt::view(Sig, i, xt::all(), xt::all()) = D;
        Sig.push_back(D);
    }
    // Sig.reverse();

    // P = mtx_norm_oracle(Sig, Y, a)
    Arr a = xt::zeros<double>({m});
    auto Q = qmi_oracle(Sig, Y);
    auto E = ell(10., a);
    auto P = bsearch_adaptor(Q, E);
    // double normY = xt::norm_l2(Y);
    auto [fb, niter, feasible] = bsearch(P, std::tuple{0., 100. * 100.});

    std::cout << niter << ", " << feasible << '\n';
    a = P.x_best();
    return niter;
    //  return prob.is_dcp()
}
