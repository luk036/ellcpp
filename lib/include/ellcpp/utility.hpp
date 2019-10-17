// -*- coding: utf-8 -*-
#pragma once

#include <type_traits>

template <typename T>
inline auto zeros(const T& x)
{
    if constexpr (std::is_floating_point_v<T>)
        return T(0.);
    else
        return T {xt::zeros<double>({x.size()})};
}
