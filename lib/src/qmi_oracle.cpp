#include <ellcpp/oracles/qmi_oracle.hpp>
#include <xtensor-blas/xlinalg.hpp>

#define myview(X, idx) xt::view(X, idx, xt::all())

using Arr = xt::xarray<double>;

/**
 * @brief
 *
 * @param x
 * @return std::tuple<Arr, double, bool>
 */
auto qmi_oracle::operator()(const Arr &x) -> std::tuple<Arr, double, bool> {
    using xt::linalg::dot;
    // using xt::placeholders::_;

    this->_count = 0;
    this->_nx = x.shape()[0];

    auto getA = [&, this](std::size_t i, std::size_t j) -> double { // ???
        assert(i >= j);
        if (this->_count < i + 1) {
            this->_count = i + 1;
            myview(this->_Fx, i) = myview(this->_F0, i);
            for (auto k = 0U; k < this->_nx; ++k) {
                myview(this->_Fx, i) -= myview(this->_F[k], i) * x(k);
            }
        }
        auto a = -dot(myview(this->_Fx, i), myview(this->_Fx, j))();
        if (i == j) {
            a += this->_t;
        }
        return a;
    };

    this->_Q.factor(getA);
    auto g = Arr{xt::zeros<double>({this->_nx})};
    if (this->_Q.is_spd()) {
        return {std::move(g), -1., true};
    }

    auto [v, ep] = this->_Q.witness();
    auto p = v.size();
    Arr Fxp = myview(this->_Fx, xt::range(0, p));
    Arr Av = dot(v, Fxp);
    for (auto k = 0U; k < this->_nx; ++k) {
        Arr Fkp = myview(this->_F[k], xt::range(0, p));
        g(k) = -2 * dot(dot(v, Fkp), Av)();
    }
    return {std::move(g), ep, false};
}
