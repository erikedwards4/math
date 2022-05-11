//Sorti/Select Help function.
//Implements algorithm to select kth largest or smallest element in a vector.
//See qselect.c for usage.
//This operates in-place and partial sorts X.

#pragma once

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


size_t lomuto_partitioni_s (FLT_I *X, const size_t N, size_t p, const int largest)
{
    const size_t hi = N - 1u;
    FLT_I x = X[p];
    float pivot = x.val;
    X[p] = X[hi]; X[hi] = x;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i].val<=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i].val>=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    x = X[p]; X[p] = X[hi]; X[hi] = x;
    return p;
}


size_t lomuto_partitioni_d (DBL_I *X, const size_t N, size_t p, const int largest)
{
    const size_t hi = N - 1u;
    DBL_I x = X[p];
    double pivot = x.val;
    X[p] = X[hi]; X[hi] = x;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i].val<=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i].val>=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    x = X[p]; X[p] = X[hi]; X[hi] = x;
    return p;
}


size_t lomuto_partitioni_c (CFLT_I *X, const size_t N, size_t p, const int largest)
{
    const size_t hi = N - 1u;
    CFLT_I x = X[p];
    float pivot = x.val;
    X[p] = X[hi]; X[hi] = x;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i].val<=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i].val>=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    x = X[p]; X[p] = X[hi]; X[hi] = x;
    return p;
}


size_t lomuto_partitioni_z (CDBL_I *X, const size_t N, size_t p, const int largest)
{
    const size_t hi = N - 1u;
    CDBL_I x = X[p];
    double pivot = x.val;
    X[p] = X[hi]; X[hi] = x;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i].val<=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i].val>=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    x = X[p]; X[p] = X[hi]; X[hi] = x;
    return p;
}


#ifdef __cplusplus
}
}
#endif
