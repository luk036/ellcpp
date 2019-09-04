// -*- coding: utf-8 -*-
#pragma once

#include <limits>
#include <xtensor/xarray.hpp>

// from itertools import chain

/*!
 * @brief Oracle for lowpass filter design
 *
 */
class lowpass_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;

  private:
    mutable unsigned int _i_Anr {0};
    mutable unsigned int _i_As {0};
    mutable unsigned int _i_Ap {0};
    // mutable unsigned int _count{0};

    const Arr& _Ap;
    const Arr& _As;
    const Arr& _Anr;
    double _Lpsq;
    double _Upsq;

  public:
    /*!
     * @brief Construct a new lowpass oracle object
     *
     * @param Ap
     * @param As
     * @param Anr
     * @param Lpsq
     * @param Upsq
     */
    lowpass_oracle(
        const Arr& Ap, const Arr& As, const Arr& Anr, double Lpsq, double Upsq)
        : _Ap {Ap}
        , _As {As}
        , _Anr {Anr}
        , _Lpsq {Lpsq}
        , _Upsq {Upsq}
    {
    }

    /*!
     * @brief
     *
     * @param x
     * @param Spsq
     * @return auto
     */
    auto operator()(const Arr& x, double Spsq) const
        -> std::tuple<Arr, Arr, double>;
};
