#include "gp_solve.hpp"

using Vec = std::valarray<double>;

template <>
Info4EM<Vec> gp_base<double>::operator()(const Vec& x) const
{
    double f;
    Vec g(x.size());

    for (auto i = 1; i != _M.size(); ++i)
    {
        g = _M[i].lse_gradient(x, f);
        if (f > 0)
            return {false, g, f, x};
    }

    g = _M[0].lse_gradient(x, f);
    return {true, g, f, x};
}
