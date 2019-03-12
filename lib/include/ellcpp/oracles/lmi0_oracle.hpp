#ifndef _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI0_ORACLE_HPP
#define _HOME_UBUNTU_CUBSTORE_ELLCPP_LMI0_ORACLE_HPP 1

//#include "mat.hpp"
#include "chol_ext.hpp"
#include <xtensor/xarray.hpp>

/**
 * @brief  Oracle for Linear Matrix Inequality
 *
 * Oracle for:
 *    F * x >= 0
 */
class lmi0_oracle {
    using Arr = xt::xarray<double>;
    using shape_type = Arr::shape_type;

  private:
    const std::vector<Arr> &_F;
    size_t _n;

  public:
    chol_ext _Q;

  public:
    /**
     * @brief Construct a new lmi0 oracle object
     *
     * @param F
     */
    explicit lmi0_oracle(const std::vector<Arr> &F)
        : _F{F},               //
          _n{F[0].shape()[0]}, //
          _Q(_n)               //
    {}

    /**
     * @brief
     *
     * @param x
     * @return auto
     */
    auto operator()(const Arr &x) -> std::tuple<Arr, double, bool>;
};

#endif