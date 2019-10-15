#include <cassert>
#include <ellcpp/oracles/qmi_oracle.hpp>
#include <xtensor-blas/xlinalg.hpp>

#define myview(X, idx) xt::view(X, idx, xt::all())

using Arr = xt::xarray<double, xt::layout_type::row_major>;

/*!
 * @brief
 *
 * @param x
 * @return std::tuple<Arr, double, bool>
 */
std::tuple<bool, Arr, double> qmi_oracle::operator()(const Arr& x)
{
    using xt::linalg::dot;

    this->_count = 0;
    this->_nx = x.shape()[0];

    auto getA = [&, this](std::size_t i, std::size_t j) -> double { // ???
        assert(i >= j);
        if (this->_count < i + 1)
        {
            this->_count = i + 1;
            myview(this->_Fx, i) = myview(this->_F0, i);
            for (size_t k = 0U; k < this->_nx; ++k)
            {
                myview(this->_Fx, i) -= myview(this->_F[k], i) * x(k);
            }
        }
        auto a = -dot(myview(this->_Fx, i), myview(this->_Fx, j))();
        if (i == j)
        {
            a += this->_t;
        }
        return a;
    };

    this->_Q.factor(getA);
    if (this->_Q.is_spd())
    {
        return {false, Arr {0.}, -1.};
    }

    auto ep = this->_Q.witness();
    // auto stop = this->_Q.stop;
    // auto start = this->_Q.start;
    auto [start, stop] = this->_Q.p;
    auto v = xt::view(this->_Q.v, xt::range(start, stop));
    auto Fxp = myview(this->_Fx, xt::range(start, stop));
    auto Av = dot(v, Fxp);
    auto g = Arr {xt::zeros<double>({this->_nx})};
    for (size_t k = 0U; k < this->_nx; ++k)
    {
        auto Fkp = myview(this->_F[k], xt::range(start, stop));
        g(k) = -2 * dot(dot(v, Fkp), Av)();
    }
    return {true, std::move(g), ep};
}
