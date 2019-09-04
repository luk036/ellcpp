// -*- coding: utf-8 -*-
#pragma once

#include <cassert>
#include <cmath>
#include <numeric>
#include <type_traits>

namespace fun
{

/*!
 * @brief Greatest common divider
 *
 * @tparam _Mn
 * @param __m
 * @param __n
 * @return _Mn
 */
template <typename _Mn>
constexpr _Mn gcd(_Mn __m, _Mn __n)
{
    return __m == 0 ? abs(__n) : __n == 0 ? abs(__m) : gcd(__n, __m % __n);
}

/*!
 * @brief Least common multiple
 *
 * @tparam _Mn
 * @param __m
 * @param __n
 * @return _Mn
 */
template <typename _Mn>
constexpr _Mn lcm(_Mn __m, _Mn __n)
{
    return (__m != 0 && __n != 0) ? (abs(__m) / gcd(__m, __n)) * abs(__n) : 0;
}

/*!
 * @brief
 *
 * @tparam Z
 */
template <typename Z>
struct Fraction
{
    using _Self = Fraction<Z>;

    Z _numerator;
    Z _denominator;

  public:
    /*!
     * @brief Construct a new Fraction object
     *
     * @param numerator
     * @param denominator
     */
    constexpr Fraction(const Z& numerator, const Z& denominator)
        : _numerator {numerator}
        , _denominator {denominator}
    {
        const Z& common = gcd(numerator, denominator);
        _numerator /= common;
        _denominator /= common;
    }

    /*!
     * @brief Construct a new Fraction object
     *
     * @param numerator
     */
    explicit Fraction(const Z& numerator)
        : _numerator {numerator}
        , _denominator {1}
    {
    }

    /*!
     * @brief Construct a new Fraction object
     *
     */
    constexpr Fraction() = default;

    /*!
     * @brief
     *
     * @return const Z&
     */
    constexpr const Z& numerator() const
    {
        return _numerator;
    }

    /*!
     * @brief
     *
     * @return const Z&
     */
    constexpr const Z& denominator() const
    {
        return _denominator;
    }

    /*!
     * @brief
     *
     * @return _Self
     */
    constexpr _Self abs() const
    {
        return _Self(std::abs(_numerator), std::abs(_denominator));
    }

    /*!
     * @brief
     *
     */
    constexpr void reciprocal()
    {
        std::swap(_numerator, _denominator);
    }

    /*!
     * @brief
     *
     * @return _Self
     */
    constexpr _Self operator-() const
    {
        return _Self(-_numerator, _denominator);
    }

    /*!
     * @brief
     *
     * @param frac
     * @return _Self
     */
    constexpr _Self operator+(const _Self& frac) const
    {
        auto common = lcm(_denominator, frac._denominator);
        auto n = common / _denominator * _numerator +
            common / frac._denominator * frac._numerator;
        return _Self(n, common);
    }

    /*!
     * @brief
     *
     * @param frac
     * @return _Self
     */
    constexpr _Self operator-(const _Self& frac) const
    {
        return *this + (-frac);
    }

    /*!
     * @brief
     *
     * @param frac
     * @return _Self
     */
    constexpr _Self operator*(const _Self& frac) const
    {
        auto n = _numerator * frac._numerator;
        auto d = _denominator * frac._denominator;
        return _Self(n, d);
    }

    /*!
     * @brief
     *
     * @param frac
     * @return _Self
     */
    constexpr _Self operator/(_Self frac) const
    {
        frac.reciprocal();
        return *this * frac;
    }

    /*!
     * @brief
     *
     * @param i
     * @return _Self
     */
    constexpr _Self operator+(const Z& i) const
    {
        auto n = _numerator + _denominator * i;
        return _Self(n, _denominator);
    }

    /*!
     * @brief
     *
     * @param i
     * @return _Self
     */
    constexpr _Self operator-(const Z& i) const
    {
        return *this + (-i);
    }

    /*!
     * @brief
     *
     * @param i
     * @return _Self
     */
    constexpr _Self operator*(const Z& i) const
    {
        auto n = _numerator * i;
        return _Self(n, _denominator);
    }

    /*!
     * @brief
     *
     * @param i
     * @return _Self
     */
    constexpr _Self operator/(const Z& i) const
    {
        auto d = _denominator * i;
        return _Self(_numerator, d);
    }

    /*!
     * @brief
     *
     * @param frac
     * @return _Self
     */
    constexpr _Self operator+=(const _Self& frac)
    {
        return *this = *this + frac;
    }

    /*!
     * @brief
     *
     * @param frac
     * @return _Self
     */
    constexpr _Self operator-=(const _Self& frac)
    {
        return *this = *this - frac;
    }

    /*!
     * @brief
     *
     * @param frac
     * @return _Self
     */
    constexpr _Self operator*=(const _Self& frac)
    {
        return *this = *this * frac;
    }

    /*!
     * @brief
     *
     * @param frac
     * @return _Self
     */
    constexpr _Self operator/=(const _Self& frac)
    {
        return *this = *this / frac;
    }

    /*!
     * @brief
     *
     * @param i
     * @return _Self
     */
    constexpr _Self operator+=(const Z& i)
    {
        return *this = *this + i;
    }

    /*!
     * @brief
     *
     * @param i
     * @return _Self
     */
    constexpr _Self operator-=(const Z& i)
    {
        return *this = *this - i;
    }

    /*!
     * @brief
     *
     * @param i
     * @return _Self
     */
    constexpr _Self operator*=(const Z& i)
    {
        return *this = *this * i;
    }

    /*!
     * @brief
     *
     * @param i
     * @return _Self
     */
    constexpr _Self operator/=(const Z& i)
    {
        return *this = *this / i;
    }

    /*!
     * @brief Three way comparison
     *
     * @tparam U
     * @param frac
     * @return auto
     */
    template <typename U>
    constexpr auto cmp(const Fraction<U>& frac) const
    {
        return _numerator * frac._denominator - _denominator * frac._numerator;
    }

    /*!
     * @brief
     *
     * @tparam U
     * @param frac
     * @return true
     * @return false
     */
    template <typename U>
    constexpr bool operator==(const Fraction<U>& frac) const
    {
        return this->cmp(frac) == 0;
    }

    /*!
     * @brief
     *
     * @tparam U
     * @param frac
     * @return true
     * @return false
     */
    template <typename U>
    constexpr bool operator!=(const Fraction<U>& frac) const
    {
        return this->cmp(frac) != 0;
    }

    /*!
     * @brief
     *
     * @tparam U
     * @param frac
     * @return true
     * @return false
     */
    template <typename U>
    constexpr bool operator<(const Fraction<U>& frac) const
    {
        return this->cmp(frac) < 0;
    }

    /*!
     * @brief
     *
     * @tparam U
     * @param frac
     * @return true
     * @return false
     */
    template <typename U>
    constexpr bool operator>(const Fraction<U>& frac) const
    {
        return this->cmp(frac) > 0;
    }

    /*!
     * @brief
     *
     * @tparam U
     * @param frac
     * @return true
     * @return false
     */
    template <typename U>
    constexpr bool operator<=(const Fraction<U>& frac) const
    {
        return this->cmp(frac) <= 0;
    }

    /*!
     * @brief
     *
     * @tparam U
     * @param frac
     * @return true
     * @return false
     */
    template <typename U>
    constexpr bool operator>=(const Fraction<U>& frac) const
    {
        return this->cmp(frac) >= 0;
    }

    /*!
     * @brief
     *
     * @param c
     * @return auto
     */
    constexpr auto cmp(const Z& c) const
    {
        return _numerator - _denominator * c;
    }

    /*!
     * @brief
     *
     * @param c
     * @return true
     * @return false
     */
    constexpr bool operator==(const Z& c) const
    {
        return this->cmp(c) == 0;
    }

    /*!
     * @brief
     *
     * @param c
     * @return true
     * @return false
     */
    constexpr bool operator!=(const Z& c) const
    {
        return this->cmp(c) != 0;
    }

    /*!
     * @brief
     *
     * @param c
     * @return true
     * @return false
     */
    constexpr bool operator<(const Z& c) const
    {
        return this->cmp(c) < 0;
    }

    /*!
     * @brief
     *
     * @param c
     * @return true
     * @return false
     */
    constexpr bool operator>(const Z& c) const
    {
        return this->cmp(c) > 0;
    }

    /*!
     * @brief
     *
     * @param c
     * @return true
     * @return false
     */
    constexpr bool operator<=(const Z& c) const
    {
        return this->cmp(c) <= 0;
    }

    /*!
     * @brief
     *
     * @param c
     * @return true
     * @return false
     */
    constexpr bool operator>=(const Z& c) const
    {
        return this->cmp(c) >= 0;
    }

    /*!
     * @brief
     *
     * @return double
     */
    explicit operator double()
    {
        return double(_numerator) / _denominator;
    }
};

/*!
 * @brief
 *
 * @tparam Z
 * @param c
 * @param frac
 * @return Fraction<Z>
 */
template <typename Z>
constexpr Fraction<Z> operator+(const Z& c, const Fraction<Z>& frac)
{
    return frac + c;
}

/*!
 * @brief
 *
 * @tparam Z
 * @param c
 * @param frac
 * @return Fraction<Z>
 */
template <typename Z>
constexpr Fraction<Z> operator-(const Z& c, const Fraction<Z>& frac)
{
    return c + (-frac);
}

/*!
 * @brief
 *
 * @tparam Z
 * @param c
 * @param frac
 * @return Fraction<Z>
 */
template <typename Z>
constexpr Fraction<Z> operator*(const Z& c, const Fraction<Z>& frac)
{
    return frac * c;
}

/*!
 * @brief
 *
 * @tparam Z
 * @param c
 * @param frac
 * @return Fraction<Z>
 */
template <typename Z>
constexpr Fraction<Z> operator+(int c, const Fraction<Z>& frac)
{
    return frac + c;
}

/*!
 * @brief
 *
 * @tparam Z
 * @param c
 * @param frac
 * @return Fraction<Z>
 */
template <typename Z>
constexpr Fraction<Z> operator-(int c, const Fraction<Z>& frac)
{
    return c + (-frac);
}

/*!
 * @brief
 *
 * @tparam Z
 * @param c
 * @param frac
 * @return Fraction<Z>
 */
template <typename Z>
constexpr Fraction<Z> operator*(int c, const Fraction<Z>& frac)
{
    return frac * c;
}

/*!
 * @brief
 *
 * @tparam _Stream
 * @tparam Z
 * @param os
 * @param frac
 * @return _Stream&
 */
template <typename _Stream, typename Z>
_Stream& operator<<(_Stream& os, const Fraction<Z>& frac)
{
    os << frac.numerator() << "/" << frac.denominator();
    return os;
}

// For template deduction
// Integral{Z} Fraction(const Z &, const Z &)->Fraction<Z>;

} // namespace fun
