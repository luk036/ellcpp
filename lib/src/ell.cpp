#include <cmath>
#include <ellcpp/ell.hpp>
#include <tuple>

ell::return_t ell::calc_ll_core(double b0, double b1, double tsq) const {
    auto b1sq = b1 * b1;
    if (b1sq > tsq || !this->_use_parallel_cut) {
        return this->calc_dc(b0, tsq);
    }

    auto params = std::tuple{0., 0., 0.};

    if (unlikely(b1 < b0)) {
        return {1, params}; // no sol'n
    }

    if (b0 == 0.) {
        return this->calc_ll_cc(b1, b1sq, tsq);
    } 

    auto n = this->_n;
    auto b0b1 = b0 * b1;
    if (unlikely(n * b0b1 < -tsq)) {
        return {3, params}; // no effect
    }

    auto t0 = tsq - b0 * b0;
    auto t1 = tsq - b1sq;
    auto bav = (b0 + b1) / 2;
    auto temp = n * bav * (b1 - b0);
    auto xi = std::sqrt( t0*t1 + temp*temp);
    auto sigma = ( n + (tsq-b0b1-xi)/(2*bav*bav) ) / (n+1);
    auto rho = sigma * bav;
    auto delta = _c1 * ((t0 + t1)/2 + xi/n) / tsq;
    params = {rho, sigma, delta};

    return {0, params};
}

/** Situation when feasible cut. */
ell::return_t ell::calc_ll_cc(double b1, double b1sq, double tsq) const {
    auto n = this->_n;
    auto temp = n*b1sq/2;
    auto xi = std::sqrt(tsq * (tsq - b1sq) + temp*temp);
    auto sigma = (n + 2*(tsq - xi) / b1sq) / (n + 1);
    auto rho = sigma * b1 / 2;
    auto delta = this->_c1 * (tsq - b1sq/2 - xi/n) / tsq;
    return {0, ell::params_t{rho, sigma, delta}};
}

/**
 * @brief Deep Cut
 */
ell::return_t ell::calc_dc(double beta, double tsq) const {
    auto params = std::tuple{0., 0., 0.};
    auto tau = std::sqrt(tsq);

    if (beta > tau) {
        return {1, params}; // no sol'n
    }

    if (beta == 0.) {
        return this->calc_cc(tsq);
    }

    auto n = this->_n;
    auto gamma = tau + n * beta;
    if (gamma < 0) {
        return {3, params}; // no effect
    }

    double rho = gamma / (n + 1);
    double sigma = 2. * rho / (tau + beta);
    double delta = this->_c1 * (tsq - beta*beta) / tsq;
    params = {rho, sigma, delta};
    return {0, params};
}

/**
 * @brief Central Cut
 */
ell::return_t ell::calc_cc(double tsq) const {
    auto np1 = this->_n + 1;
    auto sigma = 2. / np1;
    auto rho = std::sqrt(tsq) / np1;
    auto delta = _c1;
    return {0, ell::params_t{rho, sigma, delta}};
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
        return {0, tsq};
    }
    if (beta > tau) {
        return {1, tsq}; // no sol'n
    }
    if (beta < -tau) {
        return {3, tsq}; // no effect
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
    return {0, tsq};
}
