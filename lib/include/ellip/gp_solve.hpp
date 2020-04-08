#ifndef _GP_SOLVE_HPP
#define _GP_SOLVE_HPP

#include "ellipsoid.hpp"
#include "posynomial.hpp"
#include <valarray>
#include <vector>

template <typename _Tp>
class gp_base
{
    using Vec = std::valarray<double>;

  public:
    gp_base() { }
    ~gp_base() { }

    Info4EM<Vec> operator()(const Vec& x) const;

  protected:
    std::vector<posynomial<_Tp>> _M;
};

#endif
