#include <ellcpp/oracles/lmi_old_oracle.hpp>
// #include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

/*!
 * @brief
 *
 * @param x
 * @return std::tuple<Arr, double, bool>
 */
std::tuple<Arr, double> lmi_old_oracle::operator()(const Arr& x)
{
    auto n = x.size();

    auto A = Arr {this->_F0};
    for (size_t k = 0U; k < n; ++k)
    {
        A -= this->_F[k] * x(k);
    }

    this->_Q.factorize(A);
    if (this->_Q.is_spd())
    {
        return std::tuple {Arr {0.}, -1.};
    }
    auto ep = this->_Q.witness();
    auto g = Arr {xt::zeros<double>({n})};
    for (size_t i = 0U; i < n; ++i)
    {
        g(i) = this->_Q.sym_quad(this->_F[i]);
    }
    return std::tuple {std::move(g), ep};
}
