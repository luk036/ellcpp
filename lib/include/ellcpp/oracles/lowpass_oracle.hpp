// -*- coding: utf-8 -*-
#pragma once

#include <limits>
#include <xtensor/xarray.hpp>

// from itertools import chain

/*!
 * @brief Oracle for FIR lowpass filter design.
 *
 *    This example is taken from Almir Mutapcic in 2006:
 *
 *        min   \gamma
 *        s.t.  L^2(\omega) \le R(\omega) \le U^2(\omega), \forall \omega \in
 * [0, \pi] R(\omega) > 0, \forall \omega \in [0, \pi]
 */
class lowpass_oracle
{
    using Arr = xt::xarray<double, xt::layout_type::row_major>;
    using ParallelCut = std::tuple<Arr, Arr>;

  private:
    mutable size_t _i_Anr {};
    mutable size_t _i_As {};
    mutable size_t _i_Ap {};
    // mutable unsigned int _count{};

    const Arr& _Ap;
    const Arr& _As;
    const Arr& _Anr;
    double _Lpsq;
    double _Upsq;

  public:
    /*!
     * @brief Construct a new lowpass oracle object
     *
     * @param[in] Ap
     * @param[in] As
     * @param[in] Anr
     * @param[in] Lpsq
     * @param[in] Upsq
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
     * @param[in] x
     * @param[in] Spsq
     * @return auto
     */
    auto operator()(const Arr& x, double& Spsq) const
        -> std::tuple<ParallelCut, bool>;
};
