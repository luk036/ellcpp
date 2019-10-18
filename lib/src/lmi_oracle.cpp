#include <ellcpp/oracles/lmi_oracle.hpp>
#include <ellcpp/utility.hpp>
// #include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;
using Cut = std::tuple<Arr, double>;

/*!
 * @brief
 *
 * @param x
 * @return std::optional<Cut>
 */
std::optional<Cut> lmi_oracle::operator()(const Arr& x)
{
    auto n = x.size();

    auto getA = [&, this](unsigned i, unsigned j) -> double {
        auto a = this->_F0(i, j);
        for (size_t k = 0U; k < n; ++k)
        {
            a -= this->_F[k](i, j) * x(k);
        }
        return a;
    };

    this->_Q.factor(getA);
    if (this->_Q.is_spd())
    {
        return {};
    }
    auto ep = this->_Q.witness();
    auto g = zeros(x);
    for (size_t i = 0U; i < n; ++i)
    {
        g(i) = this->_Q.sym_quad(this->_F[i]);
    }
    return {{std::move(g), ep}};
}
