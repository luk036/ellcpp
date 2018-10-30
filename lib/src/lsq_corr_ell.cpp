// -*- coding: utf-8 -*-
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/lmi0_oracle.hpp>
#include <ellcpp/oracles/qmi_oracle.hpp>
#include <iostream>
#include <tuple>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xrandom.hpp>
#include <vector>

using Arr = xt::xarray<double>;


/**
 * @brief 
 * 
 * @param s 
 * @return Arr 
 */
std::tuple<Arr, Arr> create_2d_isotropic(size_t nx = 10u, size_t ny = 8u, size_t N = 3000u) {
    using xt::linalg::dot;

    const auto n = nx * ny;
    const double s_end[] = {10., 8.};
    const auto sdkern = 0.3;  // width of kernel
    const auto var = 2.;      // standard derivation
    const auto tau = 0.00001; // standard derivation of white noise
    xt::random::seed(5);

    // create sites s
    Arr sx = xt::linspace<double>(0., 10., nx);
    Arr sy = xt::linspace<double>(0., 8., ny);
    auto [xx, yy] = xt::meshgrid(sx, sy);
    Arr s = xt::stack(xt::xtuple(xt::flatten(xx), xt::flatten(yy)), 0);
    s = xt::transpose(s);

    Arr Sig = xt::zeros<double>({n, n});
    for (auto i = 0u; i < n; ++i) {
        for (auto j = i; j < n; ++j) {
            Arr d = xt::view(s, j, xt::all()) - xt::view(s, i, xt::all());
            double g = -sdkern * std::sqrt(dot(d, d)());
            Sig(i, j) = std::exp(g);
            Sig(j, i) = Sig(i, j);
        }
    }

    Arr A = xt::linalg::cholesky(Sig);
    Arr Y = xt::zeros<double>({n, n});
    for (auto k = 0u; k < N; ++k) {
        Arr x = var * xt::random::randn<double>({n});
        // auto y = dot(A, x)() + tau*xt::random::randn<double>({n});
        Arr y = dot(A, x);
        Y += xt::linalg::outer(y, y);
    }
    Y /= N;

    return std::tuple{std::move(Y), std::move(s)};
}

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
    return std::move(D);
}

class lsq_oracle {
    using Arr = xt::xarray<double>;
    using shape_type = Arr::shape_type;

  private:
    qmi_oracle _qmi;
    lmi0_oracle _lmi0;

  public:
    explicit lsq_oracle(const std::vector<Arr> &F, const Arr &F0)
        : _qmi(F, F0), //
          _lmi0(F) {}

    auto operator()(const Arr &x, double t) {
        auto n = x.shape()[0];
        Arr g = xt::zeros<double>({n});

        auto [g10, fj0, feasible0] = this->_lmi0(xt::view(x, xt::range(0, n-1)));
        if (!feasible0) {
            xt::view(g, xt::range(0, n-1)) = g10;
            g(n-1) = 0.;
            return std::tuple{std::move(g), fj0, t};
        }

        this->_qmi.update(x(n-1));
        auto [g1, fj, feasible] = this->_qmi(xt::view(x, xt::range(0, n-1)));
        if (!feasible) {
            xt::view(g, xt::range(0, n-1)) = g1;
            auto v = this->_qmi._Q.witness();
            g(n-1) = -xt::linalg::dot(v,v)();
            return std::tuple{std::move(g), fj, t};
        }

        g(n-1) = 1.;
        auto tc = x(n-1);
        auto f0 = tc - t;
        if (f0 > 0) {
            return std::tuple{std::move(g), f0, t};
        }

        return std::tuple{std::move(g), 0., tc};
    }
};


auto lsq_corr_core2(const Arr &Y, std::size_t m, lsq_oracle &P) {
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
    return std::tuple{num_iters, std::move(a)};
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

    auto P = lsq_oracle(Sig, Y);
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
    auto I = std::tuple{0., 100. * 100.};
    auto [fb, niter, feasible] = bsearch(P, I);

    std::cout << niter << ", " << feasible << '\n';
    a = P.x_best();
    return niter;
    //  return prob.is_dcp()
}
