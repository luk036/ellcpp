#include "ellipsoid.hpp"

void ellipsoid::update(const Vec& g)
{
    Vec gt(_n);
    for (size_t i = 0; i < _n; ++i)
    {
        gt[i] = (_Ae[i] * g).sum();
    }
    double gamma = (g * gt).sum();
    assert(gamma >= 0);

    double tau = sqrt(gamma);
    double sigma = 2.0 / (_n + 1) / gamma;
    double n2 = double(_n) * _n;
    double delta = n2 / (n2 - 1);

    _x -= gt / ((_n + 1) * tau);

    for (size_t i = 0; i < _n; ++i)
    {
        const double temp = sigma * gt[i];
        for (size_t j = i; j < _n; ++j)
        {
            _Ae[i][j] -= temp * gt[j];
            _Ae[i][j] *= delta;
        }
    }

    // Make symmetric
    for (size_t i = 0; i < _n - 1; ++i)
    {
        for (size_t j = i + 1; j < _n; ++j)
        {
            _Ae[j][i] = _Ae[i][j];
        }
    }
}


CUTSTATUS ellipsoid::update(const Vec& g, double beta)
{
    Vec gt(_n);
    for (size_t i = 0; i < _n; ++i)
    {
        gt[i] = (_Ae[i] * g).sum();
    }
    // b += (g*_x0).sum();
    double gamma = (g * gt).sum();
    assert(gamma >= 0);

    double tau = sqrt(gamma);
    if (beta > tau)
        return NOSOLUTION;
    if (beta < -tau)
        return NOEFFECT;

    double rho = (tau + _n * beta) / (_n + 1);
    double sigma = 2 * rho / (tau + beta) / gamma;
    double n2 = double(_n) * _n;
    double delta = n2 / (n2 - 1) * (gamma - beta * beta) / gamma;
    // double k = sigma/gamma;

    _x -= (rho / gamma) * gt;

    for (size_t i = 0; i < _n; ++i)
    {
        const double temp = sigma * gt[i];
        for (size_t j = i; j < _n; ++j)
        {
            _Ae[i][j] -= temp * gt[j];
            _Ae[i][j] *= delta;
        }
    }

    // Make symmetric
    for (size_t i = 0; i < _n - 1; ++i)
    {
        for (size_t j = i + 1; j < _n; ++j)
        {
            _Ae[j][i] = _Ae[i][j];
        }
    }

    return CUT;
}
