#ifndef _HOME_UBUNTU_GITHUB_PGCPP_FRACTION_HPP
#define _HOME_UBUNTU_GITHUB_PGCPP_FRACTION_HPP 1

#include <cassert>
#include <cmath>
#include <numeric>
#include <type_traits>

namespace fun {

template <typename _Mn> constexpr _Mn gcd(_Mn __m, _Mn __n) {
    return __m == 0 ? abs(__n) : __n == 0 ? abs(__m) : gcd(__n, __m % __n);
}

/// Least common multiple
template <typename _Mn> constexpr _Mn lcm(_Mn __m, _Mn __n) {
    return (__m != 0 && __n != 0) ? (abs(__m) / gcd(__m, __n)) * abs(__n) : 0;
}

template <typename Z>
struct Fraction {
    using _Self = Fraction<Z>;

    Z _numerator;
    Z _denominator;

  public:
    constexpr Fraction(const Z &numerator, const Z &denominator) {
        const Z &common = gcd(numerator, denominator);
        _numerator = numerator / common;
        _denominator = denominator / common;
    }

    constexpr explicit Fraction(const Z &numerator)
        : _numerator{numerator}, _denominator{1} {}

    constexpr Fraction() = default;

    constexpr const Z &numerator() const { return _numerator; }

    constexpr const Z &denominator() const { return _denominator; }

    constexpr _Self abs() const {
        return _Self(std::abs(_numerator), std::abs(_denominator));
    }

    constexpr void reciprocal() {
        std::swap(_numerator, _denominator);
    }

    constexpr _Self operator-() const {
        return _Self(-_numerator, _denominator);
    }

    constexpr _Self operator+(const _Self &frac) const {
        auto common = lcm(_denominator, frac._denominator);
        auto n = common / _denominator * _numerator +
                 common / frac._denominator * frac._numerator;
        return _Self(n, common);
    }

    constexpr _Self operator-(const _Self &frac) const {
        return *this + (-frac);
    }

    constexpr _Self operator*(const _Self &frac) const {
        auto n = _numerator * frac._numerator;
        auto d = _denominator * frac._denominator;
        return _Self(n, d);
    }

    constexpr _Self operator/(_Self frac) const {
        frac.reciprocal();
        return *this * frac;
    }

    constexpr _Self operator+(const Z &i) const {
        auto n = _numerator + _denominator * i;
        return _Self(n, _denominator);
    }

    constexpr _Self operator-(const Z &i) const { return *this + (-i); }

    constexpr _Self operator*(const Z &i) const {
        auto n = _numerator * i;
        return _Self(n, _denominator);
    }

    constexpr _Self operator/(const Z &i) const {
        auto d = _denominator * i;
        return _Self(_numerator, d);
    }

    constexpr _Self operator+=(const _Self &frac) {
        return *this = *this + frac;
    }

    constexpr _Self operator-=(const _Self &frac) {
        return *this = *this - frac;
    }

    constexpr _Self operator*=(const _Self &frac) {
        return *this = *this * frac;
    }

    constexpr _Self operator/=(const _Self &frac) {
        return *this = *this / frac;
    }

    constexpr _Self operator+=(const Z &i) { return *this = *this + i; }

    constexpr _Self operator-=(const Z &i) { return *this = *this - i; }

    constexpr _Self operator*=(const Z &i) { return *this = *this * i; }

    constexpr _Self operator/=(const Z &i) { return *this = *this / i; }

    /**
     * @brief Three way comparison
     *
     * @param frac
     * @return constexpr auto
     */
    template <typename U> constexpr auto cmp(const Fraction<U> &frac) const {
        return _numerator * frac._denominator - _denominator * frac._numerator;
    }

    template <typename U>
    constexpr bool operator==(const Fraction<U> &frac) const {
        return this->cmp(frac) == 0;
    }

    template <typename U>
    constexpr bool operator!=(const Fraction<U> &frac) const {
        return this->cmp(frac) != 0;
    }

    template <typename U>
    constexpr bool operator<(const Fraction<U> &frac) const {
        return this->cmp(frac) < 0;
    }

    template <typename U>
    constexpr bool operator>(const Fraction<U> &frac) const {
        return this->cmp(frac) > 0;
    }

    template <typename U>
    constexpr bool operator<=(const Fraction<U> &frac) const {
        return this->cmp(frac) <= 0;
    }

    template <typename U>
    constexpr bool operator>=(const Fraction<U> &frac) const {
        return this->cmp(frac) >= 0;
    }

    constexpr auto cmp(const Z &c) const {
        return _numerator - _denominator * c;
    }

    constexpr bool operator==(const Z &c) const { return this->cmp(c) == 0; }

    constexpr bool operator!=(const Z &c) const { return this->cmp(c) != 0; }

    constexpr bool operator<(const Z &c) const { return this->cmp(c) < 0; }

    constexpr bool operator>(const Z &c) const { return this->cmp(c) > 0; }

    constexpr bool operator<=(const Z &c) const { return this->cmp(c) <= 0; }

    constexpr bool operator>=(const Z &c) const { return this->cmp(c) >= 0; }

    explicit operator double() { return double(_numerator) / _denominator; }
};

template <typename Z>
constexpr Fraction<Z> operator+(const Z &c, const Fraction<Z> &frac) {
    return frac + c;
}

template <typename Z>
constexpr Fraction<Z> operator-(const Z &c, const Fraction<Z> &frac) {
    return c + (-frac);
}

template <typename Z>
constexpr Fraction<Z> operator*(const Z &c, const Fraction<Z> &frac) {
    return frac * c;
}

template <typename Z>
constexpr Fraction<Z> operator+(int c, const Fraction<Z> &frac) {
    return frac + c;
}

template <typename Z>
constexpr Fraction<Z> operator-(int c, const Fraction<Z> &frac) {
    return c + (-frac);
}

template <typename Z>
constexpr Fraction<Z> operator*(int c, const Fraction<Z> &frac) {
    return frac * c;
}

template <typename _Stream, typename Z>
_Stream &operator<<(_Stream &os, const Fraction<Z> &frac) {
    os << frac.numerator() << "/" << frac.denominator();
    return os;
}

// For template deduction
// Integral{Z} Fraction(const Z &, const Z &)->Fraction<Z>;

} // namespace fun

#endif
