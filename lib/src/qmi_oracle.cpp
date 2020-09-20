#include <cassert>
#include <ellcpp/oracles/qmi_oracle.hpp>
#include <xtensor-blas/xlinalg.hpp>

#define ROW(X, index) xt::view(X, index, xt::all())
#define COLUMN(X, index) xt::view(X, xt::all(), index)

using Arr = xt::xarray<double, xt::layout_type::row_major>;
using Cut = std::tuple<Arr, double>;

/*!
 * @brief
 *
 * @param[in] x
 * @return boost::optional<Cut>
 */
boost::optional<Cut> qmi_oracle::operator()(const Arr& x)
{
    using xt::linalg::dot;

    this->_count = 0;
    this->_nx = x.shape()[0];

    auto getA = [&, this](size_t i, size_t j) -> double { // ???
        assert(i >= j);
        if (this->_count < i + 1)
        {
            this->_count = i + 1;
            ROW(this->_Fx, i) = COLUMN(this->_F0, i);
            for (size_t k = 0U; k != this->_nx; ++k)
            {
                ROW(this->_Fx, i) -= COLUMN(this->_F[k], i) * x(k);
            }
        }
        auto a = -dot(ROW(this->_Fx, i), ROW(this->_Fx, j))();
        if (i == j)
        {
            a += this->_t;
        }
        return a;
    };

    this->_Q.factor(getA);
    if (this->_Q.is_spd())
    {
        return {};
    }

    const auto ep = this->_Q.witness();
    const auto [start, stop] = this->_Q.p;
    const auto v = xt::view(this->_Q.v, xt::range(start, stop));
    const auto Fxp = ROW(this->_Fx, xt::range(start, stop));
    const auto Av = dot(v, Fxp);
    auto g = zeros({this->_nx});
    for (size_t k = 0U; k != this->_nx; ++k)
    {
        const auto Fkp = ROW(this->_F[k], xt::range(start, stop));
        g(k) = -2 * dot(dot(v, Fkp), Av)();
    }
    return {{std::move(g), ep}};
}
