#include <ellcpp/oracles/profit_oracle.hpp>
#include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double>;

/**
 * @brief
 *
 * @param y
 * @param t
 * @return std::tuple<Arr, double, double>
 */
auto profit_oracle::operator()(const Arr &y, double t) const
    -> std::tuple<Arr, double, double> {
    // auto fj = y[0] - this->_log_k; // constraint
    if (auto fj = y[0] - this->_log_k; fj > 0.) {
        auto g = Arr{1., 0.};
        return {std::move(g), fj, t};
    }

    auto log_Cobb = this->_log_pA + xt::linalg::dot(this->_a, y)();
    auto x = Arr{xt::exp(y)};
    auto vx = xt::linalg::dot(this->_v, x)();
    auto te = t + vx;

    auto fj = std::log(te) - log_Cobb;
    if (fj < 0.) {
        te = std::exp(log_Cobb);
        t = te - vx;
        fj = 0.;
    }
    auto g = Arr{(this->_v * x) / te - this->_a};
    return {std::move(g), fj, t};
}
