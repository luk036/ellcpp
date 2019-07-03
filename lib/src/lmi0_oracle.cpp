#include <ellcpp/oracles/lmi0_oracle.hpp>
// #include <xtensor-blas/xlinalg.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

/*!
 * @brief
 *
 * @param x
 * @return auto
 */
auto lmi0_oracle::operator()(const Arr& x) -> std::tuple<Arr, double, bool>
{
    auto n = x.size();

    auto getA = [&, this](unsigned i, unsigned j) -> double {
        auto a = 0.;
        for (auto k = 0U; k < n; ++k)
        {
            a += this->_F[k](i, j) * x(k);
        }
        return a;
    };

    auto g = Arr{xt::zeros<double>({n})};

    this->_Q.factor(getA);
    if (this->_Q.is_spd()) { return std::tuple{std::move(g), -1., true}; }
    auto ep = this->_Q.witness();
    for (auto i = 0U; i < n; ++i)
    {
        g(i) = -_Q.sym_quad(this->_F[i]);
    }
    return std::tuple{std::move(g), ep, false};
}
