#include <ellcpp/oracles/profit_oracle.hpp>
#include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;
using Cut = std::tuple<Arr, double>;

/*!
 * @brief
 *
 * @param[in] y
 * @param[in] t the best-so-far optimal value
 * @return std::tuple<Cut, double>
 */
std::tuple<Cut, bool> profit_oracle::operator()(const Arr& y, double& t) const
{
    // y0 <= log k
    const auto f1 = y[0] - this->_log_k;
    if (f1 > 0.)
    {
        return {{Arr {1., 0.}, f1}, false};
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
        auto g = Arr {(this->_v * x) / te - this->_a};
        return {{std::move(g), 0.}, true};
    }
    auto g = Arr {(this->_v * x) / te - this->_a};
    return {{std::move(g), fj}, false};
}

/*!
 * @param[in] y
 * @param[in] t the best-so-far optimal value
 * @return std::tuple<Cut, double, Arr, int>
 */
std::tuple<Cut, Arr, bool, bool> profit_q_oracle::operator()(
    const Arr& y, double& t, bool retry)
{
    if (!retry)
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
        this->_yd = xt::log(x);
    }
    auto [cut, shrunk] = this->_P(this->_yd, t);
    auto& [g, h] = cut;
    h += xt::linalg::dot(g, this->_yd - y)();
    return {std::move(cut), this->_yd, shrunk, !retry};
}
