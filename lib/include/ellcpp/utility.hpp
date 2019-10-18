// -*- coding: utf-8 -*-
#pragma once

#include <type_traits>
#include <xtensor/xarray.hpp>

using Arr = xt::xarray<double, xt::layout_type::row_major>;

template <typename T>
inline auto zeros(std::initializer_list<T>&& x) -> Arr
{
    return Arr {xt::zeros<double>(x)};
}

template <typename T>
inline auto zeros(const T& x) -> T
{
    if constexpr (std::is_floating_point_v<T>)
    {
        return T(0);
    }
    else
    {
        return T {xt::zeros<double>({x.size()})};
    }
}
