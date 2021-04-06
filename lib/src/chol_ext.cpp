// -*- coding: utf-8 -*-
#include <ellcpp/oracles/chol_ext.hpp>
#include <stdexcept>


/*!
 * @brief witness that certifies $A$ is not
 * symmetric positive definite (spd)
 *
 * @return auto
 */
auto chol_ext::witness() -> double
{
    if (this->is_spd())
    {
        XTENSOR_THROW(std::runtime_error, "Implementation Error.");
    }
    // const auto& [start, n] = this->p;
    const auto& start = this->p.first;
    const auto& n = this->p.second;
    auto m = n - 1; // assume stop > 0
    this->v(m) = 1.;
    for (auto i = m; i > start; --i)
    {
        auto s = 0.;
        for (auto k = i; k <= m; ++k)
        {
            s += this->T(k, i - 1) * this->v(k);
        }
        this->v(i - 1) = -s;
    }
    return -this->T(m, m);
}

/*!
 * @brief Calculate v'*{A}(p,p)*v
 *
 * @param[in] A
 * @return double
 */
auto chol_ext::sym_quad(const Vec& A) const -> double
{
    auto res = double {};
    const auto& v = this->v;
    // const auto& [start, stop] = this->p;
    const auto& start = this->p.first;
    const auto& stop = this->p.second;
    for (auto i = start; i != stop; ++i)
    {
        auto s = double {};
        for (auto j = i + 1; j != stop; ++j)
        {
            s += A(i, j) * v(j);
        }
        res += v(i) * (A(i, i) * v(i) + 2 * s);
    }
    return res;
}

/**
 * @brief Return upper triangular matrix $R$ where $A = R^T R$
 * 
 * @return Mat 
 */
auto chol_ext::sqrt() -> Mat
{
    if (!this->is_spd())
    {
        XTENSOR_THROW(std::runtime_error, "Implementation Error.");
    }
    auto M = zeros({this->n, this->n});
    for (auto i = 0U; i != this->n; ++i)
    {
        M(i, i) = std::sqrt(this->T(i, i));
        for (auto j = i + 1; j != this->n; ++j)
        {
            M(i, j) = this->T(j, i) * M(i, i);
        }
    }
    return M;
}
