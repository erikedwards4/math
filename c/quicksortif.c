//Sortif_Help function.
//Implements quicksort algorithm (see qsortif.c for use).

#pragma once

#include "codee_math.h"
#include "insertion_sortif.c"
#include "hoare_partitionif.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


#define P_THRESH 128u


void quicksortif_s (FLT_F *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        // rearrange elements across pivot
        size_t p = hoare_partitionif_s(X, N, ascend);

        // recur on subarray containing elements that are < pivot
        if (p>P_THRESH) { quicksortif_s(X, p, ascend); }
        else { insertion_sortif_s(X, p, ascend); }

        // recur on subarray containing elements that are > pivot
        if (N-p>P_THRESH) { quicksortif_s(&X[p], N-p, ascend); }
        else { insertion_sortif_s(&X[p], N-p, ascend); }
    }
}


void quicksortif_d (DBL_D *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partitionif_d(X, N, ascend);

        if (p>P_THRESH) { quicksortif_d(X, p, ascend); }
        else { insertion_sortif_d(X, p, ascend); }

        if (N-p>P_THRESH) { quicksortif_d(&X[p], N-p, ascend); }
        else { insertion_sortif_d(&X[p], N-p, ascend); }
    }
}


void quicksortif_c (CFLT_F *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partitionif_c(X, N, ascend);

        if (p>P_THRESH) { quicksortif_c(X, p, ascend); }
        else { insertion_sortif_c(X, p, ascend); }

        if (N-p>P_THRESH) { quicksortif_c(&X[p], N-p, ascend); }
        else { insertion_sortif_c(&X[p], N-p, ascend); }
    }
}


void quicksortif_z (CDBL_D *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partitionif_z(X, N, ascend);

        if (p>P_THRESH) { quicksortif_z(X, p, ascend); }
        else { insertion_sortif_z(X, p, ascend); }

        if (N-p>P_THRESH) { quicksortif_z(&X[p], N-p, ascend); }
        else { insertion_sortif_z(&X[p], N-p, ascend); }
    }
}


#ifdef __cplusplus
}
}
#endif
