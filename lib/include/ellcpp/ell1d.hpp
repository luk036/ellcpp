// -*- coding: utf-8 -*-
#pragma once

#include <cmath>
#include <tuple>

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
        : _r {(u - l) / 2}
        , _xc {l + _r}
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
    auto xc() const -> double
    {
        return _xc;
    }

    /*!
     * @brief Set the xc object
     *
     * @param xc
     */
    void set_xc(double xc)
    {
        _xc = xc;
    }

    /*!
     * @brief
     *
     * @param g
     * @param beta
     * @return return_t
     */
    return_t update(const std::tuple<double, double>& cut);
}; // } ell1d
