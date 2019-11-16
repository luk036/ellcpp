#include <ellcpp/oracles/profit_oracle.hpp>
#include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;
using Cut = std::tuple<Arr, double>;

/*!
 * @brief
 *
 * @param y
 * @param t the best-so-far optimal value
 * @return std::tuple<Cut, double>
 */
std::tuple<Cut, double> profit_oracle::operator()(const Arr& y, double t) const
{
    // y0 <= log k
    const auto f1 = y[0] - this->_log_k;
    if (f1 > 0.)
    {
        return {{Arr {1., 0.}, f1}, t};
    }

    const auto log_Cobb = this->_log_pA + xt::linalg::dot(this->_a, y)();
    const auto x = Arr {xt::exp(y)};
    const auto vx = xt::linalg::dot(this->_v, x)();
    auto te = t + vx;

    auto fj = std::log(te) - log_Cobb;
    if (fj < 0.)
    {
        te = std::exp(log_Cobb);
        t = te - vx;
        fj = 0.;
    }
    auto g = Arr {(this->_v * x) / te - this->_a};
    return {{std::move(g), fj}, t};
}

/*!
 * @param y
 * @param t the best-so-far optimal value
 * @return std::tuple<Cut, double, Arr, int>
 */
std::tuple<Cut, double, Arr, int> profit_q_oracle::operator()(
    const Arr& y, double t, int /*unused*/) const
{
    auto x = Arr {xt::round(xt::exp(y))};
    if (x[0] == 0.)
    {
        x[0] = 1.; // nearest integer than 0
    }
    if (x[1] == 0.)
    {
        x[1] = 1.;
    }
    auto yd = Arr {xt::log(x)};
    auto [cut, t1] = this->P(yd, t);
    auto& [g, h] = cut;
    h += xt::linalg::dot(g, yd - y)();
    return {std::move(cut), t1, std::move(yd), 1};
}
