#include <ellcpp/oracles/lmi0_oracle.hpp>
// #include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

/*!
 * @brief
 *
 * @param x
 * @return auto
 */
std::tuple<Arr, double> lmi0_oracle::operator()(const Arr& x)
{
    auto n = x.size();

    auto getA = [&, this](unsigned i, unsigned j) -> double {
        auto a = 0.;
        for (size_t k = 0U; k < n; ++k)
        {
            a += this->_F[k](i, j) * x(k);
        }
        return a;
    };

    this->_Q.factor(getA);
    if (this->_Q.is_spd())
    {
        return std::tuple {Arr {0.}, -1.};
    }
    auto ep = this->_Q.witness();
    auto g = Arr {xt::zeros<double>({n})};
    for (size_t i = 0U; i < n; ++i)
    {
        g(i) = -_Q.sym_quad(this->_F[i]);
    }
    return std::tuple {std::move(g), ep};
}
