//Sortif_Help function.
//Implements quicksort algorithm (see qsortif.c for use).

#include "codee_math.h"
#include "insertion_sortif.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


#define HI_THRESH 128u

static size_t hoare_partition_s (FLT_F *X, const size_t N, const int ascend);
static size_t hoare_partition_d (DBL_D *X, const size_t N, const int ascend);
static size_t hoare_partition_c (CFLT_F *X, const size_t N, const int ascend);
static size_t hoare_partition_z (CDBL_D *X, const size_t N, const int ascend);


static size_t hoare_partition_s (FLT_F *X, const size_t N, const int ascend)
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


static size_t hoare_partition_d (DBL_D *X, const size_t N, const int ascend)
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


static size_t hoare_partition_c (CFLT_F *X, const size_t N, const int ascend)
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


static size_t hoare_partition_z (CDBL_D *X, const size_t N, const int ascend)
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


void quicksortif_s (FLT_F *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        // rearrange elements across pivot
        size_t p = hoare_partition_s(X, N, ascend);

        // recur on subarray containing elements that are < pivot
        if (p>HI_THRESH) { quicksortif_s(X, p, ascend); }
        else { insertion_sortif_s(X, p, ascend); }

        // recur on subarray containing elements that are > pivot
        if (N-p>HI_THRESH) { quicksortif_s(&X[p], N-p, ascend); }
        else { insertion_sortif_s(&X[p], N-p, ascend); }
    }
}


void quicksortif_d (DBL_D *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partition_d(X, N, ascend);

        if (p>HI_THRESH) { quicksortif_d(X, p, ascend); }
        else { insertion_sortif_d(X, p, ascend); }

        if (N-p>HI_THRESH) { quicksortif_d(&X[p], N-p, ascend); }
        else { insertion_sortif_d(&X[p], N-p, ascend); }
    }
}


void quicksortif_c (CFLT_F *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partition_c(X, N, ascend);

        if (p>HI_THRESH) { quicksortif_c(X, p, ascend); }
        else { insertion_sortif_c(X, p, ascend); }

        if (N-p>HI_THRESH) { quicksortif_c(&X[p], N-p, ascend); }
        else { insertion_sortif_c(&X[p], N-p, ascend); }
    }
}


void quicksortif_z (CDBL_D *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partition_z(X, N, ascend);

        if (p>HI_THRESH) { quicksortif_z(X, p, ascend); }
        else { insertion_sortif_z(X, p, ascend); }

        if (N-p>HI_THRESH) { quicksortif_z(&X[p], N-p, ascend); }
        else { insertion_sortif_z(&X[p], N-p, ascend); }
    }
}


#ifdef __cplusplus
}
}
#endif
