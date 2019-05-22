// -*- coding: utf-8 -*-
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/lmi0_oracle.hpp>
#include <ellcpp/oracles/lmi_oracle.hpp>
#include <ellcpp/oracles/qmi_oracle.hpp>
// #include <iostream>
#include <limits>
#include <tuple>
#include <vector>
#include <xtensor-blas/xlinalg.hpp>
// #include <xtensor/xarray.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xnorm.hpp>
#include <xtensor/xrandom.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

/**
 * @brief Create a 2d isotropic example
 *
 * @param nx
 * @param ny
 * @param N
 * @return std::tuple<Arr, Arr>
 */
std::tuple<Arr, Arr> create_2d_isotropic(size_t nx = 10u, size_t ny = 8u,
                                         size_t N = 3000u) {
    using xt::linalg::dot;

    auto const n = nx * ny;
    auto const s_end = Arr{10., 8.};
    auto const sdkern = 0.3;  // width of kernel
    auto const var = 2.;      // standard derivation
    auto const tau = 0.00001; // standard derivation of white noise
    xt::random::seed(5);

    // create sites s
    auto sx = xt::linspace<double>(0., s_end[0], nx);
    auto sy = xt::linspace<double>(0., s_end[1], ny);
    auto [xx, yy] = xt::meshgrid(sx, sy);
    auto st = Arr{xt::stack(xt::xtuple(xt::flatten(xx), xt::flatten(yy)), 0)};
    auto s = xt::transpose(st);

    auto Sig = Arr{xt::zeros<double>({n, n})};
    for (auto i = 0U; i < n; ++i) {
        for (auto j = i; j < n; ++j) {
            auto d = xt::view(s, j, xt::all()) - xt::view(s, i, xt::all());
            auto g = -sdkern * std::sqrt(dot(d, d)());
            Sig(i, j) = std::exp(g);
            Sig(j, i) = Sig(i, j);
        }
    }

    auto A = xt::linalg::cholesky(Sig);
    auto Y = Arr{xt::zeros<double>({n, n})};
    for (auto k = 0U; k < N; ++k) {
        auto x = var * xt::random::randn<double>({n});
        auto y = dot(A, x) + tau * xt::random::randn<double>({n});
        // Arr y = dot(A, x);
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
std::vector<Arr> construct_distance_matrix(const Arr &s, size_t m) {
    auto n = s.shape()[0];
    // c = cvx.Variable(m)
    auto D1 = Arr{xt::zeros<double>({n, n})};
    for (auto i = 0U; i < n; ++i) {
        for (auto j = i + 1; j < n; ++j) {
            auto h = xt::view(s, j, xt::all()) - xt::view(s, i, xt::all());
            auto d = std::sqrt(xt::linalg::dot(h, h)());
            D1(i, j) = d;
            D1(j, i) = d;
        }
    }

    auto D = Arr{xt::ones<double>({n, n})};
    auto Sig = std::vector{ D };
    Sig.reserve(m);

    for (auto i = 0U; i < m - 1; ++i) {
        D *= D1;
        Sig.push_back(D);
    }
    return std::move(Sig);
}

/**
 * @brief
 *
 */
class lsq_oracle {
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using shape_type = Arr::shape_type;

  private:
    qmi_oracle _qmi;
    lmi0_oracle _lmi0;

  public:
    /**
     * @brief Construct a new lsq oracle object
     *
     * @param F
     * @param F0
     */
    lsq_oracle(const std::vector<Arr> &F, const Arr &F0)
        : _qmi(F, F0), //
          _lmi0(F) {}

    /**
     * @brief
     *
     * @param x
     * @param t
     * @return auto
     */
    auto operator()(const Arr &x, double t) {
        auto n = x.shape()[0];
        auto g = Arr{xt::zeros<double>({n})};
        auto [g0, fj0, feasible0] =
            this->_lmi0(xt::view(x, xt::range(0, n - 1)));
        if (!feasible0) {
            xt::view(g, xt::range(0, n - 1)) = g0;
            g(n - 1) = 0.;
            return std::tuple{std::move(g), fj0, t};
        }
        this->_qmi.update(x(n - 1));

        auto [g1, fj1, feasible] = this->_qmi(xt::view(x, xt::range(0, n - 1)));
        if (!feasible) {
            xt::view(g, xt::range(0, n - 1)) = g1;
            auto &Q = this->_qmi._Q;
            auto ep = Q.witness();
            auto v = xt::view(Q.v, xt::range(Q.start, Q.stop));
            g(n - 1) = -xt::linalg::dot(v, v)();
            return std::tuple{std::move(g), fj1, t};
        }
        g(n - 1) = 1.;

        auto fj = x(n - 1) - t;
        if (fj > 0) {
            return std::tuple{std::move(g), fj, t};
        }
        return std::tuple{std::move(g), 0., x(n - 1)};
    }
};

/**
 * @brief
 *
 * @param Y
 * @param m
 * @param P
 * @return auto
 */
auto lsq_corr_core2(const Arr &Y, std::size_t m, lsq_oracle &P) {
    auto normY = 100. * xt::linalg::norm(Y);
    auto normY2 = 32. * normY * normY;
    auto val = Arr{256. * xt::ones<double>({m + 1})};

    val(m) = normY2 * normY2;
    Arr x = xt::zeros<double>({m + 1});
    x(0) = 4;
    x(m) = normY2 / 2.;
    auto E = ell(val, x);
    auto [x_best, fb, num_iters, feasible, status] =
        cutting_plane_dc(P, E, std::numeric_limits<double>::max());
    Arr a = xt::view(x_best, xt::range(0, m));
    return std::tuple{std::move(a), num_iters, feasible};
}

/**
 * @brief
 *
 * @param Y
 * @param s
 * @param m
 * @return std::tuple<size_t, bool>
 */
std::tuple<size_t, bool> lsq_corr_poly2(const Arr &Y, const Arr &s,
                                        std::size_t m) {
    auto Sig = construct_distance_matrix(s, m);
    auto P = lsq_oracle(Sig, Y);
    auto [a, num_iters, feasible] = lsq_corr_core2(Y, m, P);
    // std::cout << "lsq_corr_poly2 = " << a << "\n";
    return {num_iters, feasible};
}

/**
 * @brief
 *
 */
class mle_oracle {
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using shape_type = Arr::shape_type;

  private:
    const Arr &_Y;
    const std::vector<Arr> &_Sig;
    lmi0_oracle _lmi0;
    lmi_oracle _lmi;

  public:
    /**
     * @brief Construct a new mle oracle object
     *
     * @param Sig
     * @param Y
     */
    mle_oracle(const std::vector<Arr> &Sig, const Arr &Y)
        : _Y{Y},      //
          _Sig{Sig},  //
          _lmi0(Sig), //
          _lmi(Sig, std::move(2 * Y)) {}

    /**
     * @brief
     *
     * @param x
     * @param t
     * @return auto
     */
    auto operator()(const Arr &x, double t) {
        using xt::linalg::dot;

        auto [g1, fj1, feasible1] = this->_lmi(x);
        if (!feasible1) {
            return std::tuple{std::move(g1), fj1, t};
        }

        auto [g0, fj0, feasible0] = this->_lmi0(x);
        if (!feasible0) {
            return std::tuple{std::move(g0), fj0, t};
        }

        auto n = x.shape()[0];
        auto m = this->_Y.shape()[0];

        auto const &R = this->_lmi0._Q.sqrt();
        auto invR = Arr{xt::linalg::inv(R)};
        auto S = Arr{dot(invR, xt::transpose(invR))};
        auto SY = Arr{dot(S, this->_Y)};

        auto diag = xt::diagonal(R);
        auto f1 =
            double{2. * xt::sum(xt::log(diag))() + xt::linalg::trace(SY)()};
        // auto f1 = 0.;

        auto f = f1 - t;
        if (f < 0) {
            t = f1;
            f = 0.;
        }

        Arr g = xt::zeros<double>({n});

        for (auto i = 0U; i < n; ++i) {
            auto SFsi = dot(S, this->_Sig[i]);
            g(i) = xt::linalg::trace(SFsi)();
            for (auto k = 0U; k < m; ++k) {
                g(i) -= dot(xt::view(SFsi, k, xt::all()),
                            xt::view(SY, xt::all(), k))();
            }
        }
        return std::tuple{std::move(g), f, t};
    }
};

/**
 * @brief
 *
 * @param Y
 * @param m
 * @param P
 * @return auto
 */
auto mle_corr_core(const Arr &Y, std::size_t m, mle_oracle &P) {
    auto x = Arr{xt::zeros<double>({m})};
    x(0) = 4.;
    auto E = ell(500., x);
    auto [x_best, fb, num_iters, feasible, status] =
        cutting_plane_dc(P, E, std::numeric_limits<double>::max());
    return std::tuple{std::move(x_best), num_iters, feasible};
}

/**
 * @brief
 *
 * @param Y
 * @param s
 * @param m
 * @return std::tuple<size_t, bool>
 */
std::tuple<size_t, bool> mle_corr_poly(const Arr &Y, const Arr &s,
                                       std::size_t m) {
    auto Sig = construct_distance_matrix(s, m);
    auto P = mle_oracle(Sig, Y);
    auto [a, num_iters, feasible] = mle_corr_core(Y, m, P);
    // std::cout << "mle_corr_poly = " << a << "\n";
    return {num_iters, feasible};
}

/**
 * @brief
 *
 * @param Y
 * @param s
 * @param m
 * @return std::tuple<size_t, bool>
 */
std::tuple<size_t, bool> lsq_corr_poly(const Arr &Y, const Arr &s, size_t m) {
    auto Sig = construct_distance_matrix(s, m);
    // P = mtx_norm_oracle(Sig, Y, a)
    auto a = Arr{xt::zeros<double>({m})};
    auto Q = qmi_oracle(Sig, Y);
    auto E = ell(10., a);
    auto P = bsearch_adaptor(Q, E);
    // double normY = xt::norm_l2(Y);
    auto I = std::tuple{0., 100. * 100.};
    auto bs_info = bsearch(P, I);

    // std::cout << niter << ", " << feasible << '\n';
    a = P.x_best();
    return {bs_info.num_iters, bs_info.feasible};
    //  return prob.is_dcp()
}
