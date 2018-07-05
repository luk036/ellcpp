// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_LOWPASS_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_LOWPASS_ORACLE_HPP 1

#include <limits>
#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

// from itertools import chain

class lowpass_oracle {
    using Arr = xt::xarray<double>;

  private:
    const Arr &_Ap;
    const Arr &_As;
    const Arr &_Anr;
    double _Lpsq;
    double _Upsq;
    unsigned int _i_Anr;
    unsigned int _i_As;
    unsigned int _i_Ap;
    unsigned int _count;

  public:
    lowpass_oracle(const Arr &Ap, const Arr &As, const Arr &Anr, double Lpsq,
                   double Upsq)
        : _Ap{Ap}, _As{As}, _Anr{Anr}, _Lpsq{Lpsq}, _Upsq{Upsq}, _i_Anr{0},
          _i_As{0}, _i_Ap{0}, _count{0} {}

    auto operator()(const Arr &x, double Spsq) {
        using xt::linalg::dot;

        // 1. nonnegative-real constraint
        auto n = x.shape()[0];

        // case 1,
        if (x[0] < 0.) {
            Arr g = xt::zeros<double>({n});
            g[0] = -1.;
            Arr f{-x[0]};
            return std::tuple{g, f, Spsq};
        }

        // case 2,
        // 2. passband constraints
        auto N = _Ap.shape()[0];
        // for (k in chain(range(i_As, N), range(i_As))) {
        for (auto i = 0u, k = _i_Ap; i < N; ++i, ++k) {
            if (k == N) {
                k = 0; // round robin
            }
            auto v = dot(xt::view(_Ap, k, xt::all()), x)();
            if (v > _Upsq) {
                // f = v - Upsq;
                Arr g = xt::view(_Ap, k, xt::all());
                Arr f{v - _Upsq, v - _Lpsq};
                _i_Ap = k + 1;
                return std::tuple{g, f, Spsq};
            }

            if (v < _Lpsq) {
                // f = Lpsq - v;
                Arr g = -xt::view(_Ap, k, xt::all());
                Arr f{-v + _Lpsq, -v + _Upsq};
                _i_Ap = k + 1;
                return std::tuple{g, f, Spsq};
            }
        }

        // case 3,
        // 3. stopband constraint
        N = _As.shape()[0];
        // Arr w = xt::zeros<double>({N});
        auto fmax = std::numeric_limits<double>::min();
        auto imax = 0u;
        // for (k in chain(range(i_As, N), range(i_As))) {
        for (auto i = 0u, k = _i_As; i < N; ++i, ++k) {
            if (k == N)
                k = 0; // round robin
            auto v = dot(xt::view(_As, k, xt::all()), x)();
            if (v > Spsq) {
                // f = v - Spsq;
                Arr g = xt::view(_As, k, xt::all());
                // f = (v - Spsq, v);
                Arr f{v - Spsq, v};
                _i_As = k + 1; // k or k+1
                return std::tuple{g, f, Spsq};
            }

            if (v < 0) {
                // f = v - Spsq;
                Arr g = -xt::view(_As, k, xt::all());
                Arr f{-v, -v + Spsq};
                _i_As = k + 1;
                return std::tuple{g, f, Spsq};
            }

            if (v > fmax) {
                fmax = v;
                imax = k;
            }
        }

        // case 4,
        // 1. nonnegative-real constraint
        N = _Anr.shape()[0];
        for (auto k = 0u; k < N; ++k) {
            auto v = dot(xt::view(_Anr, k, xt::all()), x)();
            if (v < 0) {
                Arr f{-v};
                Arr g = -xt::view(_Anr, k, xt::all());
                // _i_Anr = k
                return std::tuple{g, f, Spsq};
            }
        }

        // Begin objective function
        // Spsq, imax = w.max(), w.argmax(); // update best so far Spsq
        Spsq = fmax;
        Arr f{0., fmax}; // ???
        // f = 0
        Arr g = xt::view(_As, imax, xt::all());
        return std::tuple{g, f, Spsq};
    }
};

#endif