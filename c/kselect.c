//Sort/Select Help function.
//Implements algorithm to select kth largest or smallest element in a vector.
//See qselect.c for usage.
//This operates in-place and partitions X.
//Profile note: size_t p = k is faster for median and quartiles,
//              although p = hi is a typical general choice.

#pragma once

#include "codee_math.h"
#include "lomuto_partition.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


float kselect_s (float *X, size_t N, size_t k, const int largest)
{
    while (N>1u)
    {
        size_t p = k;
        p = lomuto_partition_s(X, N, p, largest);
        if (k==p) { return *(X+k); }
        else if (k<p) { N = p; }
        else { ++p; X += p; N -= p; k -= p; }
    }
    return *X;
}


double kselect_d (double *X, size_t N, size_t k, const int largest)
{
    while (N>1u)
    {
        size_t p = k;
        p = lomuto_partition_d(X, N, p, largest);
        if (k==p) { return *(X+k); }
        else if (k<p) { N = p; }
        else { ++p; X += p; N -= p; k -= p; }
    }
    return *X;
}


size_t kselect_c (float *X, size_t N, size_t k, const int largest)
{
    size_t cnt = 0u;
    while (N>1u)
    {
        size_t p = k;
        p = lomuto_partition_c(X, N, p, largest);
        if (k==p) { return cnt+2u*k; }
        else if (k<p) { N = p; }
        else
        {
            ++p; X += 2u*p; cnt += 2u*p;
            N -= p; k -= p;
        }
    }
    return cnt;
}


size_t kselect_z (double *X, size_t N, size_t k, const int largest)
{
    size_t cnt = 0u;
    while (N>1u)
    {
        size_t p = k;
        p = lomuto_partition_z(X, N, p, largest);
        if (k==p) { return cnt+2u*k; }
        else if (k<p) { N = p; }
        else
        {
            ++p; X += 2u*p; cnt += 2u*p;
            N -= p; k -= p;
        }
    }
    return cnt;
}


#ifdef __cplusplus
}
}
#endif
