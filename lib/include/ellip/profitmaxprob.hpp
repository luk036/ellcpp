#ifndef _PROFIT_MAX_HPP
#define _PROFIX_MAX_HPP
#include "gp_solve.hpp"

/**
 * Profit Maximization Problem, taken from:
 *
 * [1] Aliabadi, Hossein, and Maziar Salahi. "Robust Geometric Programming
 * Approach to Profit Maximization with Interval Uncertainty." Computer
 * Science Journal of Moldova 21.1 (2013): 61.
 */
template <typename _Tp>
class profit_max : public gp_base<_Tp>
{
  public:
    profit_max();
    ~profit_max();

    /// Problem formulation
    void gp_setup()
    {
        using posy = posynomial<_Tp>;
        using mono = monomial<_Tp>;
        auto& M = gp_base<_Tp>::_M;

        // x0^-1
        M.emplace_back(mono(_Tp(1), {_Tp(-1), _Tp(0), _Tp(0)}));
        // p1 = x0 + v1*x1 + v2*x2
        auto p1 = mono(_Tp(1.0), {_Tp(1), _Tp(0), _Tp(0)}) +
            mono(v1, {_Tp(0), _Tp(1), _Tp(0)}) +
            mono(v2, {_Tp(0), _Tp(0), _Tp(1)});
        // k^-1 * x1
        M.emplace_back(mono(1.0 / k, {_Tp(0), _Tp(1), _Tp(0)}));
        // A * x1^alpha * x2^beta
        cd = mono(_Tp(A), {_Tp(0), alpha, beta});
        M.emplace_back(p1 / (p * cd));
    }

    /// objective function we want to analyse
    _Tp obj(const std::valarray<double>& z)
    {
        return p * exp(cd.lse(z)) - v1 * exp(z[1]) - v2 * exp(z[2]);
    }

  private:
    _Tp p;            /// @c market price per unit
    _Tp alpha;        /// @c output elasticity
    _Tp beta;         /// @c output elasticity
    _Tp v1;           /// @c first input price
    _Tp v2;           /// @c second input price
    _Tp k;            /// @c given constant that restricts the quantity of x1
    double A;         /// @c the scale of production
    monomial<_Tp> cd; /// @c Cobb-Douglas production function
};

#endif
