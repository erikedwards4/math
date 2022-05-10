//Sortif_Help function.
//Paritioning part of quicksortif algorithm.
//See quicksortif.c, etc. for usage.
//This operates in-place.

#pragma once

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


size_t hoare_partitionif_s (FLT_F *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    FLT_F x = X[0];
    float pivot = x.val;
    if (ascend)
    {
        while (1)
        {
            while (X[i].val<pivot) { ++i; }
            while (X[j].val>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i].val>pivot) { ++i; }
            while (X[j].val<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


size_t hoare_partitionif_d (DBL_D *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    DBL_D x = X[0];
    double pivot = x.val;
    if (ascend)
    {
        while (1)
        {
            while (X[i].val<pivot) { ++i; }
            while (X[j].val>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i].val>pivot) { ++i; }
            while (X[j].val<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


size_t hoare_partitionif_c (CFLT_F *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    CFLT_F x = X[0];
    float pivot = x.val;
    if (ascend)
    {
        while (1)
        {
            while (X[i].val<pivot) { ++i; }
            while (X[j].val>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i].val>pivot) { ++i; }
            while (X[j].val<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


size_t hoare_partitionif_z (CDBL_D *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    CDBL_D x = X[0];
    double pivot = x.val;
    if (ascend)
    {
        while (1)
        {
            while (X[i].val<pivot) { ++i; }
            while (X[j].val>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i].val>pivot) { ++i; }
            while (X[j].val<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


#ifdef __cplusplus
}
}
#endif
