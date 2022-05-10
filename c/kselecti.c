//Sorti/Select Help function.
//Implements algorithm to select kth largest or smallest element in a vector.
//See qselect.c for usage.
//This operates in-place and partial sorts X.

#pragma once

#include "codee_math.h"
#include "lomuto_partitioni.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


FLT_I kselecti_s (FLT_I *X, size_t hi, size_t k, const int largest)
{
    while (hi>0u)
    {
        size_t p = k;
        p = lomuto_partitioni_s(X, hi, p, largest);
        if (k==p) { return X[k]; }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
    return X[0];
}


DBL_I kselecti_d (DBL_I *X, size_t hi, size_t k, const int largest)
{
    while (hi>0u)
    {
        size_t p = k;
        p = lomuto_partitioni_d(X, hi, p, largest);
        if (k==p) { return X[k]; }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
    return X[0];
}


CFLT_I kselecti_c (CFLT_I *X, size_t hi, size_t k, const int largest)
{
    while (hi>0u)
    {
        size_t p = k;
        p = lomuto_partitioni_c(X, hi, p, largest);
        if (k==p) { return X[k]; }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
    return X[0];
}


CDBL_I kselecti_z (CDBL_I *X, size_t hi, size_t k, const int largest)
{
    while (hi>0u)
    {
        size_t p = k;
        p = lomuto_partitioni_z(X, hi, p, largest);
        if (k==p) { return X[k]; }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
    return X[0];
}


#ifdef __cplusplus
}
}
#endif
