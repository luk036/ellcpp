// -*- coding: utf-8 -*-
#include <ellcpp/cutting_plane.hpp>
#include <ellcpp/ell.hpp>
#include <ellcpp/oracles/lmi0_oracle.hpp>
#include <ellcpp/oracles/lmi_oracle.hpp>
#include <ellcpp/oracles/qmi_oracle.hpp>
#include <ellcpp/utility.hpp>
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

/*!
 * @brief Create a 2d sites object
 *
 * @param nx
 * @param ny
 * @return Arr
 */
Arr create_2d_sites(size_t nx = 10U, size_t ny = 8U)
{
    // const auto n = nx * ny;
    const auto s_end = Arr {10., 8.};
    const auto sx = xt::linspace<double>(0., s_end[0], nx);
    const auto sy = xt::linspace<double>(0., s_end[1], ny);
    const auto [xx, yy] = xt::meshgrid(sx, sy);
    const auto st =
        Arr {xt::stack(xt::xtuple(xt::flatten(xx), xt::flatten(yy)), 0)};
    return xt::transpose(st);
}

/*!
 * @brief Create a 2d isotropic object
 * 
 * @param s 
 * @param N 
 * @return Arr 
 */
Arr create_2d_isotropic(const Arr& s, size_t N = 3000U)
{
    using xt::linalg::dot;

    const auto n = s.shape()[0];
    const auto sdkern = 0.3;  // width of kernel
    const auto var = 2.;      // standard derivation
    const auto tau = 0.00001; // standard derivation of white noise
    xt::random::seed(5);

    auto Sig = zeros({n, n});
    for (auto i = 0U; i < n; ++i)
    {
        for (auto j = i; j < n; ++j)
        {
            auto d = xt::view(s, j, xt::all()) - xt::view(s, i, xt::all());
            auto g = -sdkern * std::sqrt(dot(d, d)());
            Sig(i, j) = std::exp(g);
            Sig(j, i) = Sig(i, j);
        }
    }

    auto A = xt::linalg::cholesky(Sig);
    auto Y = zeros({n, n});
    for (size_t k = 0U; k < N; ++k)
    {
        auto x = var * xt::random::randn<double>({n});
        auto y = dot(A, x) + tau * xt::random::randn<double>({n});
        // Arr y = dot(A, x);
        Y += xt::linalg::outer(y, y);
    }
    Y /= N;

    return Y;
}

/*!
 * @brief 
 * 
 * @param s 
 * @param m 
 * @return std::vector<Arr> 
 */
std::vector<Arr> construct_distance_matrix(const Arr& s, size_t m)
{
    auto n = s.shape()[0];
    // c = cvx.Variable(m)
    auto D1 = zeros({n, n});
    for (size_t i = 0U; i < n; ++i)
    {
        for (size_t j = i + 1; j < n; ++j)
        {
            auto h = xt::view(s, j, xt::all()) - xt::view(s, i, xt::all());
            auto d = std::sqrt(xt::linalg::dot(h, h)());
            D1(i, j) = d;
            D1(j, i) = d;
        }
    }

    auto D = Arr {xt::ones<double>({n, n})};
    auto Sig = std::vector {D};
    Sig.reserve(m);

    for (size_t i = 0U; i < m - 1; ++i)
    {
        D *= D1;
        Sig.push_back(D);
    }
    return Sig;
}

/*!
 * @brief
 *
 */
class lsq_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using shape_type = Arr::shape_type;
    using Cut = std::tuple<Arr, double>;

  private:
    qmi_oracle _qmi;
    lmi0_oracle _lmi0;

  public:
    /*!
     * @brief Construct a new lsq oracle object
     *
     * @param F
     * @param F0
     */
    lsq_oracle(const std::vector<Arr>& F, const Arr& F0)
        : _qmi(F, F0)
        , _lmi0(F)
    {
    }

    /*!
     * @brief
     *
     * @param x
     * @param t
     * @return auto
     */
    std::tuple<Cut, double> operator()(const Arr& x, double t)
    {
        auto n = x.size();
        auto g = zeros(x);
        if (auto cut0 = this->_lmi0(xt::view(x, xt::range(0, n - 1))))
        {
            auto [g0, f0] = *cut0;
            xt::view(g, xt::range(0, n - 1)) = g0;
            g(n - 1) = 0.;
            return {{std::move(g), f0}, t};
        }
        this->_qmi.update(x(n - 1));

        if (auto cut1 = this->_qmi(xt::view(x, xt::range(0, n - 1))))
        {
            auto [g1, f1] = *cut1;
            xt::view(g, xt::range(0, n - 1)) = g1;
            auto& Q = this->_qmi._Q;
            // auto ep = Q.witness();
            auto [start, stop] = Q.p;
            auto v = xt::view(Q.v, xt::range(start, stop));
            g(n - 1) = -xt::linalg::dot(v, v)();
            return {{std::move(g), f1}, t};
        }
        g(n - 1) = 1.;

        auto f0 = x(n - 1) - t;
        if (f0 > 0)
        {
            return {{std::move(g), f0}, t};
        }
        return {{std::move(g), 0.}, x(n - 1)};
    }
};

/*!
 * @brief
 *
 * @param Y
 * @param m
 * @param P
 * @return auto
 */
auto lsq_corr_core2(const Arr& Y, size_t m, lsq_oracle& P)
{
    auto normY = 100. * xt::linalg::norm(Y);
    auto normY2 = 32. * normY * normY;
    auto val = Arr {256. * xt::ones<double>({m + 1})};

    val(m) = normY2 * normY2;
    auto x = zeros({m + 1});
    x(0) = 4;
    x(m) = normY2 / 2.;
    auto E = ell(val, x);
    auto [x_best, ell_info] =
        cutting_plane_dc(P, E, std::numeric_limits<double>::max());
    Arr a = xt::view(x_best, xt::range(0, m));
    return std::tuple {std::move(a), ell_info.num_iters, ell_info.feasible};
}

/*!
 * @brief
 *
 * @param Y
 * @param s
 * @param m
 * @return std::tuple<size_t, bool>
 */
std::tuple<size_t, bool> lsq_corr_poly2(const Arr& Y, const Arr& s, size_t m)
{
    auto Sig = construct_distance_matrix(s, m);
    auto P = lsq_oracle(Sig, Y);
    auto [a, num_iters, feasible] = lsq_corr_core2(Y, m, P);
    // std::cout << "lsq_corr_poly2 = " << a << "\n";
    return {num_iters, feasible};
}

/*!
 * @brief
 *
 */
class mle_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using shape_type = Arr::shape_type;
    using Cut = std::tuple<Arr, double>;

  private:
    const Arr& _Y;
    const std::vector<Arr>& _Sig;
    lmi0_oracle _lmi0;
    lmi_oracle _lmi;

  public:
    /*!
     * @brief Construct a new mle oracle object
     *
     * @param Sig
     * @param Y
     */
    mle_oracle(const std::vector<Arr>& Sig, const Arr& Y)
        : _Y {Y}
        , _Sig {Sig}
        , _lmi0(Sig)
        , _lmi(Sig, 2 * Y)
    {
    }

    /*!
     * @brief
     *
     * @param x
     * @param t
     * @return auto
     */
    std::tuple<Cut, double> operator()(const Arr& x, double t)
    {
        using xt::linalg::dot;

        if (auto cut1 = this->_lmi(x))
        {
            return {*cut1, t};
        }

        if (auto cut0 = this->_lmi0(x))
        {
            return {*cut0, t};
        }

        auto n = x.shape()[0];
        auto m = this->_Y.shape()[0];

        const auto& R = this->_lmi0._Q.sqrt();
        auto invR = Arr {xt::linalg::inv(R)};
        auto S = Arr {dot(invR, xt::transpose(invR))};
        auto SY = Arr {dot(S, this->_Y)};

        auto diag = xt::diagonal(R);
        auto f1 =
            double {2. * xt::sum(xt::log(diag))() + xt::linalg::trace(SY)()};
        // auto f1 = 0.;

        auto f = f1 - t;
        if (f < 0)
        {
            t = f1;
            f = 0.;
        }

        auto g = zeros(x);

        for (size_t i = 0U; i < n; ++i)
        {
            auto SFsi = dot(S, this->_Sig[i]);
            g(i) = xt::linalg::trace(SFsi)();
            for (size_t k = 0U; k < m; ++k)
            {
                g(i) -= dot(
                    xt::view(SFsi, k, xt::all()), xt::view(SY, xt::all(), k))();
            }
        }
        return {{std::move(g), f}, t};
    }
};

/*!
 * @brief
 *
 * @param Y
 * @param m
 * @param P
 * @return auto
 */
auto mle_corr_core(const Arr& /*Y*/, size_t m, mle_oracle& P)
{
    auto x = zeros({m});
    x(0) = 4.;
    auto E = ell(500., x);
    auto [x_best, ell_info] =
        cutting_plane_dc(P, E, std::numeric_limits<double>::max());
    return std::tuple {
        std::move(x_best), ell_info.num_iters, ell_info.feasible};
}

/*!
 * @brief
 *
 * @param Y
 * @param s
 * @param m
 * @return std::tuple<size_t, bool>
 */
std::tuple<size_t, bool> mle_corr_poly(const Arr& Y, const Arr& s, size_t m)
{
    auto Sig = construct_distance_matrix(s, m);
    auto P = mle_oracle(Sig, Y);
    auto [a, num_iters, feasible] = mle_corr_core(Y, m, P);
    // std::cout << "mle_corr_poly = " << a << "\n";
    return {num_iters, feasible};
}

/*!
 * @brief
 *
 * @param Y
 * @param s
 * @param m
 * @return std::tuple<size_t, bool>
 */
std::tuple<size_t, bool> lsq_corr_poly(const Arr& Y, const Arr& s, size_t m)
{
    auto Sig = construct_distance_matrix(s, m);
    // P = mtx_norm_oracle(Sig, Y, a)
    auto a = zeros({m});
    auto Q = qmi_oracle(Sig, Y);
    auto E = ell(10., a);
    auto P = bsearch_adaptor(Q, E);
    // double normY = xt::norm_l2(Y);
    auto bs_info = bsearch(P, std::tuple {0., 100. * 100.});

    // std::cout << niter << ", " << feasible << '\n';
    a = P.x_best();
    return {bs_info.num_iters, bs_info.feasible};
    //  return prob.is_dcp()
}
