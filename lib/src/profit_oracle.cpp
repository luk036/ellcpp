#include <ellcpp/oracles/profit_oracle.hpp>
#include <xtensor-blas/xlinalg.hpp>

using xarray = xt::xarray<double>;

/**
 * @brief
 *
 * @param y
 * @param t
 * @return std::tuple<Arr, double, double>
 */
auto profit_oracle::operator()(const xarray &y, double t) const
    -> std::tuple<xarray, double, double> {
    // auto fj = y[0] - this->_log_k; // constraint
    if (auto fj = y[0] - this->_log_k; fj > 0.) {
        auto g = xarray{1., 0.};
        return {std::move(g), fj, t};
    }

    auto log_Cobb = this->_log_pA + xt::linalg::dot(this->_a, y)();
    auto x = xarray{xt::exp(y)};
    auto vx = xt::linalg::dot(this->_v, x)();
    auto te = t + vx;

    auto fj = std::log(te) - log_Cobb;
    if (fj < 0.) {
        te = std::exp(log_Cobb);
        t = te - vx;
        fj = 0.;
    }
    auto g = xarray{(this->_v * x) / te - this->_a};
    return {std::move(g), fj, t};
}
