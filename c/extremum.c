//Sort/Select Help function.
//Gets extremum of vector X with length N.
//This is the min if min is True, else the max.
//For complex, returns the amin or amax.
//This returns the value and is for internal C use only (like all the Help functions).

#pragma once

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


float extremum_s (const float *X, const size_t N, const int min)
{
    float m = *X++;
    if (min)
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (*X<m) { m = *X; } }
    }
    else
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (*X>m) { m = *X; } }
    }
    return m;
}


double extremum_d (const double *X, const size_t N, const int min)
{
    double m = *X++;
    if (min)
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (*X<m) { m = *X; } }
    }
    else
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (*X>m) { m = *X; } }
    }
    return m;
}


float extremum_c (const float *X, const size_t N, const int min)
{
    float xr = *X++, xi = *X++, xa = xr*xr + xi*xi;
    float mr = xr, mi = xi, ma = xa;
    if (min)
    {
        for (size_t n=1u; n<N; ++n)
        {
            xr = *X++; xi = *X++;
            xa = xr*xr + xi*xi;
            if (xa<ma) { mr = xr; mi = xi; ma = xa; }
        }
    }
    else
    {
        for (size_t n=1u; n<N; ++n)
        {
            xr = *X++; xi = *X++;
            xa = xr*xr + xi*xi;
            if (xa>ma) { mr = xr; mi = xi; ma = xa; }
        }
    }
    return ma;
}


double extremum_z (const double *X, const size_t N, const int min)
{
    double xr = *X++, xi = *X++;
    double mr = xr, mi = xi, xa = xr*xr+xi*xi, ma = xa;
    if (min)
    {
        for (size_t n=1u; n<N; ++n)
        {
            xr = *X++; xi = *X++;
            xa = xr*xr + xi*xi;
            if (xa<ma) { mr = xr; mi = xi; ma = xa; }
        }
    }
    else
    {
        for (size_t n=1u; n<N; ++n)
        {
            xr = *X++; xi = *X++;
            xa = xr*xr + xi*xi;
            if (xa>ma) { mr = xr; mi = xi; ma = xa; }
        }
    }
    return ma;
}


#ifdef __cplusplus
}
}
#endif
