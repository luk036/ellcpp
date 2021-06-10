#ifndef _ELLIPSOID_HPP
#define _ELLIPSOID_HPP

#include <cassert>
#include <valarray>

/** @return return status */
enum CUTSTATUS
{
    CUT,
    NOSOLUTION,
    NOEFFECT
};

/** Enclosing ellipsoid */
class ellipsoid
{
    using Vec = std::valarray<double>;
    using Matrix = std::valarray<Vec>;

  public:
    /** Constructor */
    ellipsoid(Vec& x, double rho)
        : _n(x.size())
        , _Ae(Vec(0., _n), _n)
        , _x(x)
    {
        assert(rho > 0);
        for (auto i = 0; i != _n; ++i)
        {
            _Ae[i][i] = rho; // initial radius
        }
    }

    /** Constructor */
    ellipsoid(Vec& x, const Vec& r)
        : _n(x.size())
        , _Ae(Vec(0., _n), _n)
        , _x(x)
    {
        for (auto i = 0; i != _n; ++i)
        {
            assert(r[i] > 0);
            _Ae[i][i] = r[i]; // initial radius
        }
    }

    ~ellipsoid() { }

    Vec& x()
    {
        return _x;
    }

    void update(const Vec& g);

    CUTSTATUS update(const Vec& g, double beta);

  private:
    size_t _n;
    Matrix _Ae;
    Vec _x; //< centroid of ellipsoid
};


/** @return return status */
enum STATUS
{
    FOUND,
    EXCEEDMAXITER,
    NOTFOUND
};


/**
 * -- (Generalized) bisection method for solving convex minimization problem P:
 *
 *     minimize     fct_0(x)
 *     subject to   fct_j(x) <= 0
 *       where fct_0(x) and fct_j(x)'s are convex
 *
 * Input:
 * 	E(x)          initial enclosing region
 *      max_it        maximum number of iterations
 *      tol           error tolerance
 *      P             Representation of convex minimization problem
 *
 * Requirement of P:
 *	void		P.assess(x)         assessment of x
 *	bool		P.is_feasible()     return true if x is feasible
 *	double		P.f_value()         returns fct_j(x) if x is infeasible for some
 *j fct::Vec	P.subgradient()     returns subgradient of fct_0 if x is
 *feasible subgradient of fct_j if x is infeasible for some j Requirement of E:
 *
 * output
 *      x             optimal solution
 *      status        FOUND = solution found to tolerance
 *                    EXCEEDMAXITER = no convergence given max_it
 *                    NOTFOUND = no feasible sol'n
 */
//#include <iostream>

template <class Vec>
struct Info4EM
{
    bool _is_feasible; /// if x is feasible
    Vec _g;            /// subgradient at x
    double _f;         /// function value at x
    Vec _x;            /// the nearest solution
};

#include <limits>
template <class Vec>
inline double norm(const Vec& x)
{
    return sqrt((x * x).sum());
}

template <class Enclosing, class Oracle, class Vec>
STATUS ellipsoid_dc_discrete(Enclosing& E, Oracle& P, Vec& best_x,
    double& best_f, int max_it = 100, double tol = 1e-4)
{
    STATUS flag = NOTFOUND;
    CUTSTATUS status = CUT;

    Vec lx = E.x();
    best_f = 1.e100; // std::numeric_limits<double>::max()

    for (int iter = 1; iter <= max_it; ++iter)
    {
        const Info4EM<Vec> info = P(lx);
        const Vec& x = info._x;
        const Vec& g = info._g;
        double f = info._f;
        double beta = (g * (lx - x)).sum();

        if (info._is_feasible)
        {
            flag = FOUND;
            if (f < best_f)
            {
                best_f = f;
                best_x = x;
            }
            status = E.update(g, beta);
        }
        else
        {
            status = E.update(g, f + beta);
        }

        if (status != CUT)
            return flag;
        if (norm(lx - E.x()) < tol)
            return flag;
        lx = E.x();
    }

    return EXCEEDMAXITER;
}


template <class Enclosing, class Oracle, class Vec>
STATUS ellipsoid_dc(Enclosing& E, Oracle& P, Vec& best_x, double& best_f,
    int max_it = 100, double tol = 1e-4)
{
    STATUS flag = NOTFOUND;
    CUTSTATUS status = CUT;

    Vec lx = E.x();
    best_f = 1.e100; // std::numeric_limits<double>::max()

    for (int iter = 1; iter <= max_it; ++iter)
    {
        const Info4EM<Vec> info = P(lx);
        // const Vec& x = info._x;
        const Vec& g = info._g;
        double f = info._f;
        // double beta = (g * (lx - x)).sum();

        if (info._is_feasible)
        {
            flag = FOUND;
            if (f < best_f)
            {
                best_f = f;
                best_x = x;
            }
            // status = E.update(g, beta);
            status = E.update(g);
        }
        else
        {
            // status = E.update(g, f + beta);
            status = E.update(g, f);
        }

        if (status != CUT)
            return flag;
        if (norm(lx - E.x()) < tol)
            return flag;
        lx = E.x();
    }

    return EXCEEDMAXITER;
}


template <class Enclosing, class Oracle, class Vec>
STATUS ellipsoid_algo(Enclosing& E, Oracle& P, Vec& best_x, double& best_f,
    int max_it = 100, double tol = 1e-4)
{
    STATUS flag = NOTFOUND;
    Vec lx = E.x();
    best_f = 1.e100; // std::numeric_limits<double>::max()

    for (int iter = 1; iter <= max_it; ++iter)
    {
        const Info4EM<Vec> info = P(lx);
        const Vec& g = info._g;
        double f = info._f;

        if (info._is_feasible)
        {
            flag = FOUND;
            if (f < best_f)
            {
                best_f = f;
                best_x = lx;
            }
        }
        E.update(g);
        if (norm(lx - E.x()) < tol)
            return flag;
        lx = E.x();
    }

    return EXCEEDMAXITER;
}


#endif
