//Sort/Select Help function.
//Implements algorithm to select kth largest or smallest element in a vector.
//See qselect.c for usage.
//This operates in-place and partial sorts X.
//Profile note: size_t p = k is faster for median and quartiles,
//but size_t p = hi is faster for deciles (and is typical general choice).

#pragma once

#include "codee_math.h"
#include "lomuto_partition.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


float kselect_s (float *X, size_t hi, size_t k, const int largest)
{
    while (hi>0u)
    {
        size_t p = k;
        p = lomuto_partition_s(X, hi, p, largest);
        if (k==p) { return *(X+k); }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
    return *X;
}


double kselect_d (double *X, size_t hi, size_t k, const int largest)
{
    while (hi>0u)
    {
        size_t p = k;
        p = lomuto_partition_d(X, hi, p, largest);
        if (k==p) { return *(X+k); }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
    return *X;
}


size_t kselect_c (float *X, size_t hi, size_t k, const int largest)
{
    size_t cnt = 0u;
    while (hi>0u)
    {
        size_t p = k;
        p = lomuto_partition_c(X, hi, p, largest);
        if (k==p) { return cnt+2u*k; }
        else if (k<p) { hi = p - 1u; }
        else
        {
            ++p; X += 2u*p; cnt += 2u*p;
            hi -= p; k -= p;
        }
    }
    return cnt;
}


size_t kselect_z (double *X, size_t hi, size_t k, const int largest)
{
    size_t cnt = 0u;
    while (hi>0u)
    {
        size_t p = k;
        p = lomuto_partition_z(X, hi, p, largest);
        if (k==p) { return cnt+2u*k; }
        else if (k<p) { hi = p - 1u; }
        else
        {
            ++p; X += 2u*p; cnt += 2u*p;
            hi -= p; k -= p;
        }
    }
    return cnt;
}


#ifdef __cplusplus
}
}
#endif
