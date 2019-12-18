// -*- coding: utf-8 -*-
#pragma once

#include <type_traits>
#include <utility>
#include <xtensor/xarray.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

template <typename T>
constexpr typename std::enable_if<std::is_integral<T>::value, Arr>::type  //
zeros(std::initializer_list<T>&& x)
{
    return Arr {xt::zeros<double>(x)};
}

template <typename T>
constexpr typename std::enable_if<std::is_floating_point<T>::value, T>::type  //
zeros(const T& x)
{
    return T{};
}

template <typename T>
constexpr typename std::enable_if<!std::is_floating_point<T>::value, T>::type  //
zeros(const T& x)
{
    return T {xt::zeros<double>({x.size()})};
}
