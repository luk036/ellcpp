// -*- coding: utf-8 -*-
#pragma once

#include <type_traits>
#include <utility>
#include <xtensor/xarray.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

template <typename T>
constexpr auto zeros(std::initializer_list<T>&& x) -> typename std::enable_if<std::is_integral<T>::value, Arr>::type
{
    return Arr {xt::zeros<double>(x)};
}

template <typename T>
constexpr auto zeros(const T& /* unused */) noexcept(noexcept(T{})) -> typename std::enable_if<std::is_floating_point<T>::value, T>::type 
{
    return T {};
}

template <typename T>
constexpr auto zeros(const T& x) -> typename std::enable_if<!std::is_floating_point<T>::value, T>::type
{
    return T {xt::zeros<double>({x.size()})};
}
