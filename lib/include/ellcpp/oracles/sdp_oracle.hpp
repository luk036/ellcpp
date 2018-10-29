// -*- coding: utf-8 -*-
#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_SDP_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_SDP_ORACLE_HPP 1

#include "lmi_oracle.hpp"
namespace bnu = boost::numeric::ublas;

class sdp_oracle {
    using Mat = bnu::symmetric_matrix<double, bnu::upper>;
    using Vec = bnu::vector<double>;
    using Arr = bnu::vector<Mat>;

  public:
    sdp_oracle(const Vec &c, const Arr &F) : _c{c}, _lmi{lmi_oracle(F)} {}

    auto operator()(const Vec &x, double t) const {
        auto [g, fj] = _lmi.chk_spd(x);
        if (fj > 0) {
            return std::make_tuple(g, fj, t);
        }
        auto f0 = bnu::inner_prod(_c, x);
        fj = f0 - t;
        g = _c;
        if (fj > 0) {
            return std::make_tuple(g, fj, t);
        }
        return std::make_tuple(g, 0., f0);
    }

  private:
    const Vec &_c;
    lmi_oracle _lmi;
};

#endif