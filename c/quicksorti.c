//Sorti_Help function.
//Implements quicksort algorithm (see qsorti.c for use).

#include "codee_math.h"
#include "insertion_sorti.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


#define HI_THRESH 128u

static size_t hoare_partition_s (FLT_I *X, const size_t N, const int ascend);
static size_t hoare_partition_d (DBL_I *X, const size_t N, const int ascend);
static size_t hoare_partition_c (CFLT_I *X, const size_t N, const int ascend);
static size_t hoare_partition_z (CDBL_I *X, const size_t N, const int ascend);


static size_t hoare_partition_s (FLT_I *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    FLT_I x = X[0];
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


static size_t hoare_partition_d (DBL_I *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    DBL_I x = X[0];
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


static size_t hoare_partition_c (CFLT_I *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    CFLT_I x = X[0];
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


static size_t hoare_partition_z (CDBL_I *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    CDBL_I x = X[0];
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


void quicksorti_s (FLT_I *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        // rearrange elements across pivot
        size_t p = hoare_partition_s(X, N, ascend);

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
        size_t p = hoare_partition_d(X, N, ascend);

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
        size_t p = hoare_partition_c(X, N, ascend);

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
        size_t p = hoare_partition_z(X, N, ascend);

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
