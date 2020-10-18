#include <ellcpp/oracles/lmi0_oracle.hpp>
#include <ellcpp/utility.hpp>
// #include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;
using Cut = std::tuple<Arr, double>;

/*!
 * @brief
 *
 * @param[in] x
 * @return auto
 */
std::optional<Cut> lmi0_oracle::operator()(const Arr& x)
{
    auto n = x.size();

    auto getA = [&, this](size_t i, size_t j) -> double {
        auto a = 0.;
        for (size_t k = 0U; k != n; ++k)
        {
            a += this->_F[k](i, j) * x(k);
        }
        return a;
    };

    if (this->_Q.factor(getA))
    {
        return {};
    }
    auto ep = this->_Q.witness();
    auto g = zeros(x);
    for (size_t i = 0U; i != n; ++i)
    {
        g(i) = -_Q.sym_quad(this->_F[i]);
    }
    return {{std::move(g), ep}};
}
