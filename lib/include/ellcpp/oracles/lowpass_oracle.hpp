// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_LOWPASS_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_LOWPASS_ORACLE_HPP 1

#include <limits>
#include <xtensor/xarray.hpp>

// from itertools import chain

/**
 * @brief Oracle for lowpass filter design
 *
 */
class lowpass_oracle {
    using Arr = xt::xarray<double>;

  private:
    const Arr &_Ap;
    const Arr &_As;
    const Arr &_Anr;
    double _Lpsq;
    double _Upsq;
    mutable unsigned int _i_Anr;
    mutable unsigned int _i_As;
    mutable unsigned int _i_Ap;
    mutable unsigned int _count;

  public:
    /**
     * @brief Construct a new lowpass oracle object
     *
     * @param Ap
     * @param As
     * @param Anr
     * @param Lpsq
     * @param Upsq
     */
    lowpass_oracle(const Arr &Ap, const Arr &As, const Arr &Anr, double Lpsq,
                   double Upsq)
        : _Ap{Ap}, _As{As}, _Anr{Anr}, _Lpsq{Lpsq}, _Upsq{Upsq}, _i_Anr{0},
          _i_As{0}, _i_Ap{0}, _count{0} {}

    /**
     * @brief
     *
     * @param x
     * @param Spsq
     * @return auto
     */
    auto operator()(const Arr &x, double Spsq) const
        -> std::tuple<Arr, Arr, double>;
};

#endif
