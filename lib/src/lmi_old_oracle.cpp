#include <ellcpp/oracles/lmi_old_oracle.hpp>
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
std::optional<Cut> lmi_old_oracle::operator()(const Arr& x)
{
    const auto n = x.size();

    auto A = Arr {this->_F0};
    for (size_t k = 0U; k != n; ++k)
    {
        A -= this->_F[k] * x(k);
    }

    this->_Q.factorize(A);
    if (this->_Q.is_spd())
    {
        return {};
    }
    const auto ep = this->_Q.witness();
    auto g = zeros(x);
    for (size_t i = 0U; i != n; ++i)
    {
        g(i) = this->_Q.sym_quad(this->_F[i]);
    }
    return {{std::move(g), ep}};
}
