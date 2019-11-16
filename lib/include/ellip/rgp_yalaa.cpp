#include "gp_solve.hpp"
#include <map>
#include <yalaa.hpp>

typedef yalaa::details::double_iv_t iv_t;
typedef yalaa::traits::interval_traits<iv_t> iv_traits;
typedef yalaa::aff_e_d aaf;
typedef std::map<typename aaf::error_t, int> pmap;

template <typename AF>
typename AF::base_t max(const AF& af, pmap& pol)
{
    typename AF::ac_t ac(af.ac());
    typename AF::base_t res = ac.central();
    for (auto it(ac.begin()); it != ac.end(); ++it)
    {
        if (it->dev() > 0)
        {
            pol[*it] = 1;
            res += it->dev();
        }
        else if (it->dev() < 0)
        {
            pol[*it] = -1;
            res -= it->dev();
        }
        else
        {
            pol[*it] = 0;
        }
    }
    return res;
}

template <typename AF>
typename AF::base_t eval(const AF& af, const pmap& pol)
{
    typename AF::ac_t ac(af.ac());
    typename AF::base_t res = ac.central();
    for (auto it(ac.begin()); it != ac.end(); ++it)
    {
        pmap::const_iterator pit(pol.find(*it));
        if (pit == pol.end())
            throw; // Noise symbol is not in the map
        res += it->dev() * pit->second;
    }
    return res;
}

/** Constructor (for aaf -> double) */
template <>
template <>
monomial<double>::monomial(const monomial<aaf>& mon, const pmap& polarity)
    : _b(eval(mon._b, polarity))
    , _a(mon._a.size())
{
    for (size_t i = 0; i < _a.size(); ++i)
    {
        _a[i] = eval(mon._a[i], polarity);
    }
}

/** Constructor (for AAF -> double) */
template <>
template <>
posynomial<double>::posynomial(
    const posynomial<aaf>& posyn, const pmap& polarity)
    : _M(posyn._M.size(), monomial<double>(posyn._M[0]._a.size()))
{
    for (size_t i = 0; i < _M.size(); ++i)
    {
        _M[i] = monomial<double>(posyn._M[i], polarity);
    }
}


using Vec = std::valarray<double>;

template <>
Info4EM<Vec> gp_base<aaf>::operator()(const Vec& x) const
{
    double f;
    Vec g(x.size());
    pmap polarity;

    for (size_t i = 1; i < _M.size(); ++i)
    {
        f = max(_M[i].lse(x), polarity);
        if (f > 0)
        {
            posynomial<double> P(_M[i], polarity);
            g = P.lse_gradient(x, f);
            if (f > 0)
            { // double check for robustness
                return {true, g, f, x};
            }
        }
    }

    f = max(_M[0].lse(x), polarity);
    posynomial<double> P(_M[0], polarity);
    g = P.lse_gradient(x, f);
    return {true, g, f, x};
}
