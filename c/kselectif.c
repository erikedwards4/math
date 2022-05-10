//Sortif/Select Help function.
//Implements algorithm to select kth largest or smallest element in a vector.
//See qselect.c for usage.
//This operates in-place and partial sorts X.

#pragma once

#include "codee_math.h"
#include "lomuto_partitionif.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


FLT_F kselectif_s (FLT_F *X, size_t hi, size_t k, const int largest)
{
    while (hi>0u)
    {
        size_t p = k;
        p = lomuto_partitionif_s(X, hi, p, largest);
        if (k==p) { return X[k]; }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
    return X[0];
}


DBL_D kselectif_d (DBL_D *X, size_t hi, size_t k, const int largest)
{
    while (hi>0u)
    {
        size_t p = k;
        p = lomuto_partitionif_d(X, hi, p, largest);
        if (k==p) { return X[k]; }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
    return X[0];
}


CFLT_F kselectif_c (CFLT_F *X, size_t hi, size_t k, const int largest)
{
    while (hi>0u)
    {
        size_t p = k;
        p = lomuto_partitionif_c(X, hi, p, largest);
        if (k==p) { return X[k]; }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
    return X[0];
}


CDBL_D kselectif_z (CDBL_D *X, size_t hi, size_t k, const int largest)
{
    while (hi>0u)
    {
        size_t p = k;
        p = lomuto_partitionif_z(X, hi, p, largest);
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
