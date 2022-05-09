//Sorti_Help function.
//Implements quicksort algorithm (see qsorti.c for use).

#include "codee_math.h"
#include "insertion_sorti.c"
#include "hoare_partitioni.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


#define HI_THRESH 128u


void quicksorti_s (FLT_I *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        // rearrange elements across pivot
        size_t p = hoare_partitioni_s(X, N, ascend);

        // recur on subarray containing elements that are < pivot
        if (p>HI_THRESH) { quicksorti_s(X, p, ascend); }
        else { insertion_sorti_s(X, p, ascend); }

        // recur on subarray containing elements that are > pivot
        if (N-p>HI_THRESH) { quicksorti_s(&X[p], N-p, ascend); }
        else { insertion_sorti_s(&X[p], N-p, ascend); }
    }
}


void quicksorti_d (DBL_I *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partitioni_d(X, N, ascend);

        if (p>HI_THRESH) { quicksorti_d(X, p, ascend); }
        else { insertion_sorti_d(X, p, ascend); }

        if (N-p>HI_THRESH) { quicksorti_d(&X[p], N-p, ascend); }
        else { insertion_sorti_d(&X[p], N-p, ascend); }
    }
}


void quicksorti_c (CFLT_I *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partitioni_c(X, N, ascend);

        if (p>HI_THRESH) { quicksorti_c(X, p, ascend); }
        else { insertion_sorti_c(X, p, ascend); }

        if (N-p>HI_THRESH) { quicksorti_c(&X[p], N-p, ascend); }
        else { insertion_sorti_c(&X[p], N-p, ascend); }
    }
}


void quicksorti_z (CDBL_I *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partitioni_z(X, N, ascend);

        if (p>HI_THRESH) { quicksorti_z(X, p, ascend); }
        else { insertion_sorti_z(X, p, ascend); }

        if (N-p>HI_THRESH) { quicksorti_z(&X[p], N-p, ascend); }
        else { insertion_sorti_z(&X[p], N-p, ascend); }
    }
}


#ifdef __cplusplus
}
}
#endif
