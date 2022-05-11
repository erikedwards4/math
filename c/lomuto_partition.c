//Sort/Select Help function.
//Partitioning part of quickselect algorithm.
//See kselect.c for usage.
//This operates in-place.

#pragma once

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


size_t lomuto_partition_s (float *X, const size_t N, size_t p, const int largest)
{
    const size_t hi = N - 1u;
    float x = X[p], pivot = x;
    X[p] = X[hi]; X[hi] = x;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i]<=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i]>=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    x = X[p]; X[p] = X[hi]; X[hi] = x;
    return p;
}


size_t lomuto_partition_d (double *X, const size_t N, size_t p, const int largest)
{
    const size_t hi = N - 1u;
    double x = X[p], pivot = x;
    X[p] = X[hi]; X[hi] = x;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i]<=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i]>=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    x = X[p]; X[p] = X[hi]; X[hi] = x;
    return p;
}


size_t lomuto_partition_c (float *X, const size_t N, size_t p, const int largest)
{
    const size_t hi = N - 1u;
    float xr = X[2u*p], xi = X[2u*p+1u], pivot = xr*xr + xi*xi;
    X[2u*p] = X[2u*hi]; X[2u*hi] = xr;
    X[2u*p+1u] = X[2u*hi+1u]; X[2u*hi+1u] = xi;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[2u*i]*X[2u*i]+X[2u*i+1u]*X[2u*i+1u]<=pivot+1e-5f)
            {
                xr = X[2u*i]; X[2u*i] = X[2u*p]; X[2u*p] = xr;
                xi = X[2u*i+1u]; X[2u*i+1u] = X[2u*p+1u]; X[2u*p+1u] = xi;
                ++p;
            }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[2u*i]*X[2u*i]+X[2u*i+1u]*X[2u*i+1u]>=pivot-1e-5f)
            {
                xr = X[2u*i]; X[2u*i] = X[2u*p]; X[2u*p] = xr;
                xi = X[2u*i+1u]; X[2u*i+1u] = X[2u*p+1u]; X[2u*p+1u] = xi;
                ++p;
            }
        }
    }
    xr = X[2u*p]; X[2u*p] = X[2u*hi]; X[2u*hi] = xr;
    xi = X[2u*p+1u]; X[2u*p+1u] = X[2u*hi+1u]; X[2u*hi+1u] = xi;
    return p;
}


size_t lomuto_partition_z (double *X, const size_t N, size_t p, const int largest)
{
    const size_t hi = N - 1u;
    double xr = X[2u*p], xi = X[2u*p+1u], pivot = xr*xr + xi*xi;
    X[2u*p] = X[2u*hi]; X[2u*hi] = xr;
    X[2u*p+1u] = X[2u*hi+1u]; X[2u*hi+1u] = xi;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[2u*i]*X[2u*i]+X[2u*i+1u]*X[2u*i+1u]<=pivot+1e-9)
            {
                xr = X[2u*i]; X[2u*i] = X[2u*p]; X[2u*p] = xr;
                xi = X[2u*i+1u]; X[2u*i+1u] = X[2u*p+1u]; X[2u*p+1u] = xi;
                ++p;
            }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[2u*i]*X[2u*i]+X[2u*i+1u]*X[2u*i+1u]>=pivot-1e-9)
            {
                xr = X[2u*i]; X[2u*i] = X[2u*p]; X[2u*p] = xr;
                xi = X[2u*i+1u]; X[2u*i+1u] = X[2u*p+1u]; X[2u*p+1u] = xi;
                ++p;
            }
        }
    }
    xr = X[2u*p]; X[2u*p] = X[2u*hi]; X[2u*hi] = xr;
    xi = X[2u*p+1u]; X[2u*p+1u] = X[2u*hi+1u]; X[2u*hi+1u] = xi;
    return p;
}


#ifdef __cplusplus
}
}
#endif
