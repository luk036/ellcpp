#include <ellcpp/oracles/lmi_oracle.hpp>
#include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double>;

/**
 * @brief
 *
 * @param x
 * @return std::tuple<Arr, double, bool>
 */
auto lmi_oracle::operator()(const Arr &x) -> std::tuple<Arr, double, bool> {
    using xt::linalg::dot;
    // using xt::placeholders::_;
    auto n = x.size();

    auto getA = [&, this](unsigned i, unsigned j) -> double {
        auto a = this->_F0(i, j);
        for (auto k = 0U; k < n; ++k) {
            // const auto &Fi = this->_F[k];
            a -= this->_F[k](i, j) * x(k);
        }
        return a;
    };

    auto g = Arr{xt::zeros<double>({n})};
    this->_Q.factor(getA);
    if (this->_Q.is_spd()) {
        return {std::move(g), -1., true};
    }
    auto [v, ep] = this->_Q.witness();
    for (auto i = 0U; i < n; ++i) {
        // const auto &Fi = this->_F[i];
        g(i) = this->_Q.sym_quad(v, this->_F[i]);
    }
    return {std::move(g), ep, false};
}
