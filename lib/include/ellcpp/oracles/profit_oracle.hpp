#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_PROFIT_ORACLE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_PROFIT_ORACLE_HPP 1

//#include <valarray>
#include <cmath>
#include <tuple>

#include <xtensor-blas/xlinalg.hpp>
#include <xtensor/xarray.hpp>

class profit_oracle {
    using xarray = xt::xarray<double>;

  private:
    double _log_pA;
    double _log_k;
    xarray _v;

  public:
    xarray _a;

  public:
    profit_oracle(double p, double A, double k, xarray a, xarray v)
        : _log_pA{std::log(p * A)}, //
          _log_k{std::log(k)},      //
          _v{v},                    //
          _a{a}                     //
    {}

    auto operator()(const xarray &y, double t) const {
        double fj = y[0] - _log_k; // constraint
        if (fj > 0.) {
            auto g = xarray{1., 0.};
            return std::tuple{g, fj, t};
        }
        double log_Cobb = _log_pA + xt::linalg::dot(_a, y)();
        auto x = xt::exp(y);
        double vx = xt::linalg::dot(_v, x)();
        auto te = t + vx;
        fj = std::log(te) - log_Cobb;
        if (fj < 0.) {
            te = std::exp(log_Cobb);
            t = te - vx;
            fj = 0.;
        }
        xarray g = (_v * x) / te - _a;
        return std::tuple{g, fj, t};
    }
};

class profit_rb_oracle {
    using xarray = xt::xarray<double>;

  private:
    xarray _uie;
    xarray _a;
    double _uie3;
    profit_oracle _P;

  public:
    profit_rb_oracle(double p, double A, double k, xarray a, xarray v,
                     double ui, xarray e, double e3)
        : _uie{ui * e},                             //
          _a{a},                                    //
          _uie3{ui * e3},                           //
          _P(p - _uie3, A, k - _uie3, a, v + _uie3) //
    {}

    auto operator()(const xarray &y, double t) {
        auto a_rb = _a;
        for (auto i : {0, 1}) {
            a_rb[i] += _uie[i] * (y[i] > 0. ? -1 : +1);
        }
        _P._a = a_rb;
        return _P(y, t);
    }
};

class profit_q_oracle {
    using xarray = xt::xarray<double>;

  private:
    profit_oracle _P;

  public:
    profit_q_oracle(double p, double A, double k, xarray a, xarray v)
        : _P(p, A, k, a, v) {}

    auto operator()(const xarray &y, double t, int /*unused*/) const {
        xarray x = xt::round(xt::exp(y));
        if (x[0] == 0.0) {
            x[0] = 1.0;
        }
        if (x[1] == 0.0) {
            x[1] = 1.0;
        }
        xarray yd = xt::log(x);
        auto [g, fj, t1] = _P(yd, t);
        return std::tuple{g, fj, t1, yd, 1};
    }
};

#endif
