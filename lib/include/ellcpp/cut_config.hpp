#pragma once

#include <cstddef>

enum class CUTStatus
{
    success,
    nosoln,
    smallenough,
    noeffect
};

/*!
 * @brief Options
 *
 */
struct Options
{
    unsigned int max_it = 2000; //!< maximum number of iterations
    double tol = 1e-8;          //!< error tolerance
};

/*!
 * @brief CInfo
 *
 */
struct CInfo
{
    bool feasible;
    size_t num_iters;
    CUTStatus status;
};
