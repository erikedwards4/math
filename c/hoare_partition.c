//Sort_Help function.
//Paritioning part of quicksort algorithm.
//See quicksort.c, etc. for usage.
//This operates in-place.

#pragma once

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


size_t hoare_partition_s (float *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    float x, pivot = *X;
    if (ascend)
    {
        while (1)
        {
            while (X[i]<pivot) { ++i; }
            while (X[j]>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i]>pivot) { ++i; }
            while (X[j]<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


size_t hoare_partition_d (double *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    double x, pivot = *X;
    if (ascend)
    {
        while (1)
        {
            while (X[i]<pivot) { ++i; }
            while (X[j]>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i]>pivot) { ++i; }
            while (X[j]<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


size_t hoare_partition_c (float *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = 2u*N - 2u;
    float xr, xi, pivot = X[0]*X[0] + X[1]*X[1];
    if (ascend)
    {
        while (1)
        {
            while (X[i]*X[i]+X[i+1u]*X[i+1u]<pivot-1e-5f) { i += 2u; }
            while (X[j]*X[j]+X[j+1u]*X[j+1u]>pivot+1e-5f) { j -= 2u; }
            if (j<=i) { return j/2u+1u; }
            xr = X[i]; X[i] = X[j]; X[j] = xr;
            xi = X[i+1u]; X[i+1u] = X[j+1u]; X[j+1u] = xi;
            i += 2u; j -= 2u;
        }
    }
    else
    {
        while (1)
        {
            while (X[i]*X[i]+X[i+1u]*X[i+1u]>pivot+1e-5f) { i += 2u; }
            while (X[j]*X[j]+X[j+1u]*X[j+1u]<pivot-1e-5f) { j -= 2u; }
            if (j<=i) { return j/2u+1u; }
            xr = X[i]; X[i] = X[j]; X[j] = xr;
            xi = X[i+1u]; X[i+1u] = X[j+1u]; X[j+1u] = xi;
            i += 2u; j -= 2u;
        }
    }
}


size_t hoare_partition_z (double *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = 2u*N - 2u;
    double xr, xi, pivot = X[0]*X[0] + X[1]*X[1];
    if (ascend)
    {
        while (1)
        {
            while (X[i]*X[i]+X[i+1u]*X[i+1u]<pivot-1e-9) { i += 2u; }
            while (X[j]*X[j]+X[j+1u]*X[j+1u]>pivot+1e-9) { j -= 2u; }
            if (j<=i) { return j/2u+1u; }
            xr = X[i]; X[i] = X[j]; X[j] = xr;
            xi = X[i+1u]; X[i+1u] = X[j+1u]; X[j+1u] = xi;
            i += 2u; j -= 2u;
        }
    }
    else
    {
        while (1)
        {
            while (X[i]*X[i]+X[i+1u]*X[i+1u]>pivot+1e-9) { i += 2u; }
            while (X[j]*X[j]+X[j+1u]*X[j+1u]<pivot-1e-9) { j -= 2u; }
            if (j<=i) { return j/2u+1u; }
            xr = X[i]; X[i] = X[j]; X[j] = xr;
            xi = X[i+1u]; X[i+1u] = X[j+1u]; X[j+1u] = xi;
            i += 2u; j -= 2u;
        }
    }
}


#ifdef __cplusplus
}
}
#endif
