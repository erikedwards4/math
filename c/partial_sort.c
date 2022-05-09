//Sort/Select Help function.
//Partial sort of elements in X up to the kth element.
//Based on quickselect algorithm.
//See partsort.c for usage.
//This operates in-place.
//The m2start functions move m (min or max) to X[0], and X[0] to index of m.

#include "codee_math.h"
#include "lomuto_partition.c"
#include "quicksort.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


#define K_THRESH 2u

static void m2start_s (float *X, size_t N, const int ascend);
static void m2start_d (double *X, size_t N, const int ascend);
static void m2start_c (float *X, size_t N, const int ascend);
static void m2start_z (double *X, size_t N, const int ascend);


static void m2start_s (float *X, size_t N, const int ascend)
{
    size_t i = 0u;
    float m = *X++;
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (*X<m) { m = *X; i = n; } }
    }
    else
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (*X>m) { m = *X; i = n; } }
    }
    X -= N; X[i] = X[0]; X[0] = m;
}


static void m2start_d (double *X, size_t N, const int ascend)
{
    size_t i = 0u;
    double m = *X++;
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (*X<m) { m = *X; i = n; } }
    }
    else
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (*X>m) { m = *X; i = n; } }
    }
    X -= N; X[i] = X[0]; X[0] = m;
}


static void m2start_c (float *X, size_t N, const int ascend)
{
    size_t i = 0u;
    float xr = *X++, xi = *X++;
    float mr = xr, mi = xi, xa = xr*xr+xi*xi, ma = xa;
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n)
        {
            xr = *X++; xi = *X++;
            xa = xr*xr + xi*xi;
            if (xa<ma) { mr = xr; mi = xi; ma = xa; i = n; }
        }
    }
    else
    {
        for (size_t n=1u; n<N; ++n)
        {
            xr = *X++; xi = *X++;
            xa = xr*xr + xi*xi;
            if (xa>ma) { mr = xr; mi = xi; ma = xa; i = n; }
        }
    }
    X -= 2u*N;
    X[2u*i] = X[0]; X[2u*i+1u] = X[1]; X[0] = mr; X[1] = mi;
}


static void m2start_z (double *X, size_t N, const int ascend)
{
    size_t i = 0u;
    double xr = *X++, xi = *X++;
    double mr = xr, mi = xi, xa = xr*xr+xi*xi, ma = xa;
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n)
        {
            xr = *X++; xi = *X++;
            xa = xr*xr + xi*xi;
            if (xa<ma) { mr = xr; mi = xi; ma = xa; i = n; }
        }
    }
    else
    {
        for (size_t n=1u; n<N; ++n)
        {
            xr = *X++; xi = *X++;
            xa = xr*xr + xi*xi;
            if (xa>ma) { mr = xr; mi = xi; ma = xa; i = n; }
        }
    }
    X -= 2u*N;
    X[2u*i] = X[0]; X[2u*i+1u] = X[1]; X[0] = mr; X[1] = mi;
}


void partial_sort_s (float *X, size_t N, size_t k, const int ascend)
{
    if (k<K_THRESH)
    {
        for (size_t j=0u; j<=k; ++j) { m2start_s(X+j, N-j, ascend); }
    }
    else
    {
        if (k<N-1u)
        {
            size_t p, cnt=0u;
            while (N>1u)
            {
                p = lomuto_partition_s(X, N-1u, N-1u, !ascend);
                if (k==p) { break; }
                else if (k<p) { N = p; }
                else { ++p; X += p; N -= p; k -= p; cnt += p; }
            }
            X -= cnt; k += cnt;
        }
        quicksort_s(X, k+1u, ascend);
    }
}


void partial_sort_d (double *X, size_t N, size_t k, const int ascend)
{
    if (k<K_THRESH)
    {
        for (size_t j=0u; j<=k; ++j) { m2start_d(X+j, N-j, ascend); }
    }
    else
    {
        if (k<N-1u)
        {
            size_t p, cnt=0u;
            while (N>1u)
            {
                p = lomuto_partition_d(X, N-1u, N-1u, !ascend);
                if (k==p) { break; }
                else if (k<p) { N = p; }
                else { ++p; X += p; N -= p; k -= p; cnt += p; }
            }
            X -= cnt; k += cnt;
        }
        quicksort_d(X, k+1u, ascend);
    }
}


void partial_sort_c (float *X, size_t N, size_t k, const int ascend)
{
    if (k<K_THRESH)
    {
        for (size_t j=0u; j<=k; ++j) { m2start_c(X+2u*j, N-j, ascend); }
    }
    else
    {
        if (k<N-1u)
        {
            size_t p, cnt=0u;
            while (N>1u)
            {
                p = lomuto_partition_c(X, N-1u, N-1u, !ascend);
                if (k==p) { break; }
                else if (k<p) { N = p; }
                else { ++p; X += 2u*p; cnt += 2u*p; N -= p; k -= p; }
            }
            X -= cnt; k += cnt;
        }
        quicksort_c(X, k+1u, ascend);
    }
}


void partial_sort_z (double *X, size_t N, size_t k, const int ascend)
{
    if (k<K_THRESH)
    {
        for (size_t j=0u; j<=k; ++j) { m2start_z(X+2u*j, N-j, ascend); }
    }
    else
    {
        if (k<N-1u)
        {
            size_t p, cnt=0u;
            while (N>1u)
            {
                p = lomuto_partition_z(X, N-1u, N-1u, !ascend);
                if (k==p) { break; }
                else if (k<p) { N = p; }
                else { ++p; X += 2u*p; cnt += 2u*p; N -= p; k -= p; }
            }
            X -= cnt; k += cnt;
        }
        quicksort_z(X, k+1u, ascend);
    }
}


#ifdef __cplusplus
}
}
#endif
