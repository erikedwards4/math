//Sort_Help function.
//Implements quicksort algorithm (see qsort.c for use).
//This operates in-place.

#include "codee_math.h"
#include "insertion_sort.c"
#include "hoare_partition.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


#define HI_THRESH 128u


void quicksort_s (float *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        // rearrange elements across pivot
        size_t p = hoare_partition_s(X, N, ascend);

        // recur on subarray containing elements that are < pivot
        if (p>HI_THRESH) { quicksort_s(X, p, ascend); }
        else { insertion_sort_s(X, p, ascend); }

        // recur on subarray containing elements that are > pivot
        if (N-p>HI_THRESH) { quicksort_s(X+p, N-p, ascend); }
        else { insertion_sort_s(X+p, N-p, ascend); }
    }
}


void quicksort_d (double *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partition_d(X, N, ascend);

        if (p>HI_THRESH) { quicksort_d(X, p, ascend); }
        else { insertion_sort_d(X, p, ascend); }

        if (N-p>HI_THRESH) { quicksort_d(X+p, N-p, ascend); }
        else { insertion_sort_d(X+p, N-p, ascend); }
    }
}


void quicksort_c (float *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partition_c(X, N, ascend);

        if (p>HI_THRESH) { quicksort_c(X, p, ascend); }
        else { insertion_sort_c(X, p, ascend); }

        if (N-p>HI_THRESH) { quicksort_c(X+2u*p, N-p, ascend); }
        else { insertion_sort_c(X+2u*p, N-p, ascend); }
    }
}


void quicksort_z (double *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partition_z(X, N, ascend);

        if (p>HI_THRESH) { quicksort_z(X, p, ascend); }
        else { insertion_sort_z(X, p, ascend); }

        if (N-p>HI_THRESH) { quicksort_z(X+2u*p, N-p, ascend); }
        else { insertion_sort_z(X+2u*p, N-p, ascend); }
    }
}


#ifdef __cplusplus
}
}
#endif
