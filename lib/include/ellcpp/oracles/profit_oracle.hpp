#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_PROFIT_ORACLE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_PROFIT_ORACLE_HPP 1

//#include <valarray>
#include <cmath>
#include <tuple>
#include <xtensor/xarray.hpp>

/**
 * @brief Oracle for a profit maximization problem
 *
 */
class profit_oracle {
    using xarray = xt::xarray<double>;

  private:
    double _log_pA;
    double _log_k;
    xarray _v;

  public:
    xarray _a;

  public:
    /**
     * @brief Construct a new profit oracle object
     *
     * @param p
     * @param A
     * @param k
     * @param a
     * @param v
     */
    profit_oracle(double p, double A, double k, xarray a, xarray v)
        : _log_pA{std::log(p * A)}, //
          _log_k{std::log(k)},      //
          _v{v},                    //
          _a{a}                     //
    {}

    /**
     * @brief
     *
     * @param y
     * @param t
     * @return auto
     */
    auto operator()(const xarray &y, double t) const
                    -> std::tuple<xarray, double, double>;
};

/**
 * @brief Oracle for a profit maximization problem (robust version)
 *
 */
class profit_rb_oracle {
    using xarray = xt::xarray<double>;

  private:
    xarray _uie;
    xarray _a;
    double _uie3;
    profit_oracle _P;

  public:
    /**
     * @brief Construct a new profit rb oracle object
     *
     * @param p
     * @param A
     * @param k
     * @param a
     * @param v
     * @param ui
     * @param e
     * @param e3
     */
    profit_rb_oracle(double p, double A, double k, xarray a, xarray v,
                     double ui, xarray e, double e3)
        : _uie{ui * e},                             //
          _a{a},                                    //
          _uie3{ui * e3},                           //
          _P(p - _uie3, A, k - _uie3, a, v + _uie3) //
    {}

    /**
     * @brief
     *
     * @param y
     * @param t
     * @return auto
     */
    auto operator()(const xarray &y, double t) {
        xarray a_rb = _a;
        for (auto &&i : {0, 1}) {
            a_rb[i] += _uie[i] * (y[i] > 0 ? -1 : +1);
        }
        _P._a = a_rb;
        return _P(y, t);
    }
};

/**
 * @brief Oracle for profit maximization problem (discrete version)
 *
 */
class profit_q_oracle {
    using xarray = xt::xarray<double>;

  private:
    profit_oracle _P;

  public:
    /**
     * @brief Construct a new profit q oracle object
     *
     * @param p
     * @param A
     * @param k
     * @param a
     * @param v
     */
    profit_q_oracle(double p, double A, double k, xarray a, xarray v)
        : _P(p, A, k, a, v) {}

    /**
     * @brief
     *
     * @param y
     * @param t
     * @return auto
     */
    auto operator()(const xarray &y, double t, int /*unused*/) const {
        xarray x = xt::round(xt::exp(y));
        if (x[0] == 0.) {
            x[0] = 1.;
        }
        if (x[1] == 0.) {
            x[1] = 1.;
        }
        xarray yd = xt::log(x);
        auto [g, fj, t1] = _P(yd, t);
        return std::tuple{std::move(g), fj, t1, yd, 1};
    }
};

#endif
