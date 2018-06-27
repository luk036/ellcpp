#include <cmath>
#include <ellcpp/ell.hpp>
#include <tuple>

/**
 * @brief Central Cut
 */
ell::params_t ell::calc_cc(double tsq) const {
    auto np1 = this->_n + 1;
    auto sigma = 2. / np1;
    double rho = std::sqrt(tsq) / np1;
    double delta = _c1;
    return std::tuple{rho, sigma, delta};
}

/**
 * @brief Deep Cut
 */
ell::return_t ell::calc_dc(double h0, double tsq) const {
    if (h0 == 0.) {
        return std::tuple{0, this->calc_cc(tsq)};
    }
    auto params = std::tuple{0., 0., 0.};

    auto t0 = tsq - h0*h0;
    if (t0 < 0) {
        return std::tuple{1, params}; // no sol'n
    }

    auto n = this->_n;
    auto tau = std::sqrt(tsq);
    auto gamma = tau + n * h0;    
    if (gamma < 0) {
        return std::tuple{3, params}; // no effect
    }

    double rho = gamma / (n + 1);
    double sigma = 2. * rho / (tau + h0);
    double delta = this->_c1 * t0/tsq;
    params = std::tuple{rho, sigma, delta};
    return std::tuple{0, params};
}

/** Situation when feasible cut. */
ell::params_t ell::calc_ll_cc(double h1, double t1, double tsq) const {
    auto n = this->_n;
    auto hsq1 = tsq - t1;
    auto temp = n * hsq1 / 2;
    auto xi = std::sqrt(tsq*t1 + temp*temp);
    auto sigma = (n + 2*(tsq - xi) / hsq1) / (n + 1);
    auto rho = sigma * h1 / 2;
    auto delta = this->_c1 * (tsq - hsq1/2 - xi/n) / tsq;
    return std::tuple{rho, sigma, delta};
}

ell::return_t ell::calc_ll_core(double h0, double h1, double tsq) const {
    auto t1 = tsq - h1*h1;
    if (t1 < 0 || !this->_use_parallel) {
        return this->calc_dc(h0, tsq);
    }
    auto params = std::tuple{0., 0., 0.};
    auto l = h1 - h0;
    if (l < 0) {
        return std::tuple{1, params}; // no sol'n
    }
    auto n = this->_n;
    auto p = h0*h1;
    if (n*p < -tsq) {
        return std::tuple{3, params}; // no effect
    }
    
    if (h0 == 0.) {
        params = this->calc_ll_cc(h1, t1, tsq);
    } else {
        auto t0 = tsq - h0*h0;
        auto h = (h0 + h1)/2;
        auto temp = n*h*l;
        auto xi = std::sqrt(t0*t1 + temp*temp);
        auto sigma = (n + (tsq - p - xi)/(2*h*h)) / (n + 1);
        auto rho = sigma * h;
        auto delta = _c1 * ((t0 + t1)/2 + xi/n) / tsq;
        params = std::tuple{rho, sigma, delta};
    }
    return std::tuple{0, params};
}


ell1d::return_t ell1d::update(double g, double beta) {
    auto tau = std::abs(_r * g);
    auto tsq = tau * tau;
    if (beta == 0.) {
        _r /= 2;
        if (g > 0.) {
            _xc -= _r;
        } else {
            _xc += _r;
        }
        return std::tuple{0, tsq};
    }
    if (beta > tau) {
        return std::tuple{1, tsq}; // no sol'n
    }
    if (beta < -tau) {
        return std::tuple{3, tsq}; // no effect
    }
    double l, u;
    double bound = _xc - beta / g;
    if (g > 0.) {
        u = bound;
        l = _xc - _r;
    } else {
        l = bound;
        u = _xc + _r;
    }
    _r = (u - l) / 2;
    _xc = l + _r;
    return std::tuple{0, tsq};
}
