#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_PROFIT_ORACLE_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_PROFIT_ORACLE_HPP 1

//#include <valarray>
#include <cmath>
#include <tuple>

#include <boost/numeric/ublas/vector.hpp>
namespace bnu = boost::numeric::ublas;

class profit_oracle {
  using Vec = bnu::vector<double>;
  // using Vec = std::valarray<double>;

public:
  profit_oracle(double p, double A, double alpha, double beta, double v1,
                double v2, double k)
      : _log_pA{std::log(p * A)}, _log_k{std::log(k)}, _v(2), _a(2) {
    _v[0] = v1;
    _v[1] = v2;
    _a[0] = alpha;
    _a[1] = beta;
  }

  auto operator()(const Vec &y, double t) const {
    double fj = y[0] - _log_k; // constraint
    if (fj > 0.0) {
      Vec g(2); // boost::numeric::ublas is too old!
      g[0] = 1.;
      g[1] = 0.;
      return std::tuple{g, fj, t};
    }
    auto log_Cobb = _log_pA + bnu::inner_prod(_a, y);
    // auto x = bnu::exp(y);
    Vec x(2);
    x[0] = std::exp(y[0]);
    x[1] = std::exp(y[1]);
    auto vx = bnu::inner_prod(_v, x);
    auto te = t + vx;
    fj = std::log(te) - log_Cobb;
    if (fj < 0.0) {
      te = std::exp(log_Cobb);
      t = te - vx;
      fj = 0.0;
    }
    Vec g = bnu::element_prod(_v, x) / te - _a;
    return std::tuple{g, fj, t};
  }

private:
  double _log_pA;
  double _log_k;
  Vec _v;
  Vec _a;
};

class profit_rb_oracle {
  // using namespace boost::numeric::ublas;
  using Vec = bnu::vector<double>;
  // using Vec = std::valarray<double>;

public:
  profit_rb_oracle(double p, double A, double alpha, double beta, double v1,
                   double v2, double k, double ui, double e1, double e2,
                   double e3)
      : _uie1{ui * e1}, _uie2{ui * e2}, _log_pA{std::log((p - ui * e3) * A)},
        _log_k{std::log(k - ui * e3)}, _v(2), _a(2) {
    _v[0] = v1 + ui * e3;
    _v[1] = v2 + ui * e3;
    _a[0] = alpha;
    _a[1] = beta;
  }

  auto operator()(const Vec &y, double t) const {
    double fj = y[0] - _log_k; // constraint
    if (fj > 0.0) {
      Vec g(2);
      g[0] = 1.;
      g[1] = 0.;
      return std::tuple{g, fj, t};
    }
    auto a_rb = _a;
    a_rb[0] += _uie1 * (y[0] > 0. ? -1 : +1);
    a_rb[1] += _uie2 * (y[1] > 0. ? -1 : +1);

    auto log_Cobb = _log_pA + bnu::inner_prod(a_rb, y);
    // auto x = bnu::exp(y);
    Vec x(2);
    x[0] = std::exp(y[0]);
    x[1] = std::exp(y[1]);
    auto vx = bnu::inner_prod(_v, x);
    auto te = t + vx;
    fj = std::log(te) - log_Cobb;
    if (fj < 0.0) {
      te = std::exp(log_Cobb);
      t = te - vx;
      fj = 0.0;
    }
    Vec g = bnu::element_prod(_v, x) / te - a_rb;
    return std::tuple{g, fj, t};
  }

private:
  double _uie1, _uie2;
  double _log_pA;
  double _log_k;
  Vec _v;
  Vec _a;
};

class profit_q_oracle : public profit_oracle {
  using Vec = bnu::vector<double>;
  // using Vec = std::valarray<double>;

public:
  profit_q_oracle(double p, double A, double alpha, double beta, double v1,
                  double v2, double k)
      : profit_oracle(p, A, alpha, beta, v1, v2, k) {}

  auto operator()(const Vec &y, double t, int /*unused*/) const {
    Vec yd = y;
    for (double &e : yd) {
      auto temp = std::round(std::exp(e));
      if (temp == 0.0) {
        temp = 1.0;
      }
      e = std::log(temp);
    }
    auto [g, fj, t1] = profit_oracle::operator()(yd, t);
    return std::tuple{g, fj, t1, yd, 1};
  }
};

#endif
