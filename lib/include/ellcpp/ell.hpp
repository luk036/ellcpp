#ifndef _HOME_UBUNTU_GITHUB_ELLCPP_ELL_HPP
#define _HOME_UBUNTU_GITHUB_ELLCPP_ELL_HPP 1

#include <cmath>
#include <tuple>
#include <xtensor/xarray.hpp>

/*!
 * @brief Ellipsoid Search Space
 *
 * ell = { x | (x - xc)' * (\kappa Q)^-1 * (x - xc) <= 1 }
 */
class ell
{
public:
    using Arr      = xt::xarray<double, xt::layout_type::row_major>;
    using params_t = std::tuple<double, double, double>;
    using return_t = std::tuple<int, params_t>;

public:
    bool _use_parallel_cut = true;

private:
    std::size_t _n;
    double      _c1;
    double      _kappa;
    Arr         _xc;
    Arr         _Q;

public:
    /*!
     * @brief Construct a new ell object
     *
     * @tparam T
     * @param val
     * @param x
     */
    ell(const Arr& val, const Arr& x)
        : _n{x.size()},                         // n
          _c1{double(_n * _n) / (_n * _n - 1)}, //
          _xc{x}
    {
        this->_Q     = xt::diag(val);
        this->_kappa = 1.;
    }

    /*!
     * @brief Construct a new ell object
     *
     * @tparam T
     * @param val
     * @param x
     */
    ell(const double& alpha, const Arr& x)
        : _n{x.size()},                         // n
          _c1{double(_n * _n) / (_n * _n - 1)}, //
          _xc{x}
    {
        this->_Q     = xt::eye(_n);
        this->_kappa = alpha;
    }

    /*!
     * @brief Construct a new ell object
     *
     * @param E
     */
    ell(const ell& E) = default;

    /*!
     * @brief
     *
     * @return const Arr&
     */
    auto xc() const -> const Arr& { return _xc; }

    /*!
     * @brief
     *
     * @return Arr&
     */
    auto xc() -> Arr& { return _xc; }

    /*!
     * @brief Set the xc object
     *
     * @param xc
     */
    void set_xc(const Arr& xc) { _xc = xc; }

    /*!
     * @brief Update ellipsoid core function using the cut
     *          g' * (x - xc) + beta <= 0
     *
     * @tparam T
     * @param g
     * @param beta
     * @return std::tuple<int, double>
     */
    template<typename T>
    std::tuple<int, double> update(const Arr& g, const T& beta);

private:
    /*!
     * @brief
     *
     * @param b0
     * @param b1
     * @param tsq
     * @return return_t
     */
    return_t calc_ll_core(double b0, double b1, double tsq) const;

    /*!
     * @brief Parallel Cut, one of them is central
     *
     * @param b1
     * @param t1
     * @param tsq
     * @return return_t
     */
    return_t calc_ll_cc(double b1, double t1, double tsq) const;

    /*!
     * @brief Deep Cut
     *
     * @param b0
     * @param tsq
     * @return return_t
     */
    return_t calc_dc(double b0, double tsq) const;

    /*!
     * @brief Central Cut
     *
     * @param tsq
     * @return return_t
     */
    return_t calc_cc(double tsq) const;
}; // } ell

/*!
 * @brief Ellipsoid Method for special 1D case
 *
 */
class ell1d
{
public:
    using return_t = std::tuple<int, double>;

private:
    double _r;
    double _xc;

public:
    /*!
     * @brief Construct a new ell1d object
     *
     * @param l
     * @param u
     */
    ell1d(double l, double u) //
        : _r{(u - l) / 2},    //
          _xc{l + _r}
    {
    }

    /*!
     * @brief Construct a new ell1d object
     *
     * @param E
     */
    ell1d(const ell1d& E) = default;

    /*!
     * @brief
     *
     * @return double
     */
    double xc() const { return _xc; }

    /*!
     * @brief
     *
     * @param g
     * @param beta
     * @return return_t
     */
    return_t update(double g, double beta);
}; // } ell1d

#endif
