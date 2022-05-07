//Sort_Help function.
//Implements quicksort algorithm (see qsort.c for use).

#include "codee_math.h"
#include "insertion_sort.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


#define HI_THRESH 128u

static size_t hoare_partition_s (float *X, const size_t N, const int ascend);
static size_t hoare_partition_d (double *X, const size_t N, const int ascend);
static size_t hoare_partition_c (float *X, const size_t N, const int ascend);
static size_t hoare_partition_z (double *X, const size_t N, const int ascend);


static size_t hoare_partition_s (float *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    float x, pivot = *X;
    if (ascend)
    {
        while (1)
        {
            while (X[i]<pivot) { ++i; }
            while (X[j]>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i]>pivot) { ++i; }
            while (X[j]<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


static size_t hoare_partition_d (double *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    double x, pivot = *X;
    if (ascend)
    {
        while (1)
        {
            while (X[i]<pivot) { ++i; }
            while (X[j]>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i]>pivot) { ++i; }
            while (X[j]<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


static size_t hoare_partition_c (float *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = 2u*N - 2u;
    float xr, xi, pivot = X[0]*X[0] + X[1]*X[1];
    if (ascend)
    {
        while (1)
        {
            while (X[i]*X[i]+X[i+1u]*X[i+1u]<pivot-1e-5f) { i += 2u; }
            while (X[j]*X[j]+X[j+1u]*X[j+1u]>pivot+1e-5f) { j -= 2u; }
            if (j<=i) { return j/2u+1u; }
            xr = X[i]; X[i] = X[j]; X[j] = xr;
            xi = X[i+1u]; X[i+1u] = X[j+1u]; X[j+1u] = xi;
            i += 2u; j -= 2u;
        }
    }
    else
    {
        while (1)
        {
            while (X[i]*X[i]+X[i+1u]*X[i+1u]>pivot+1e-5f) { i += 2u; }
            while (X[j]*X[j]+X[j+1u]*X[j+1u]<pivot-1e-5f) { j -= 2u; }
            if (j<=i) { return j/2u+1u; }
            xr = X[i]; X[i] = X[j]; X[j] = xr;
            xi = X[i+1u]; X[i+1u] = X[j+1u]; X[j+1u] = xi;
            i += 2u; j -= 2u;
        }
    }
}


static size_t hoare_partition_z (double *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = 2u*N - 2u;
    double xr, xi, pivot = X[0]*X[0] + X[1]*X[1];
    if (ascend)
    {
        while (1)
        {
            while (X[i]*X[i]+X[i+1u]*X[i+1u]<pivot-1e-9) { i += 2u; }
            while (X[j]*X[j]+X[j+1u]*X[j+1u]>pivot+1e-9) { j -= 2u; }
            if (j<=i) { return j/2u+1u; }
            xr = X[i]; X[i] = X[j]; X[j] = xr;
            xi = X[i+1u]; X[i+1u] = X[j+1u]; X[j+1u] = xi;
            i += 2u; j -= 2u;
        }
    }
    else
    {
        while (1)
        {
            while (X[i]*X[i]+X[i+1u]*X[i+1u]>pivot+1e-9) { i += 2u; }
            while (X[j]*X[j]+X[j+1u]*X[j+1u]<pivot-1e-9) { j -= 2u; }
            if (j<=i) { return j/2u+1u; }
            xr = X[i]; X[i] = X[j]; X[j] = xr;
            xi = X[i+1u]; X[i+1u] = X[j+1u]; X[j+1u] = xi;
            i += 2u; j -= 2u;
        }
    }
}


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
