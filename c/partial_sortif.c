//Sortif/Select Help function.
//Partial sort of elements in X up to the kth element.
//Based on quickselect algorithm.
//This operates in-place.
//The m2start functions move m (min or max) to X[0], and X[0] to index of m.

#pragma once

#include "codee_math.h"
#include "lomuto_partitionif.c"
#include "quicksortif.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


#define K_THRESH 5u

static void m2startif_s (FLT_F *X, size_t N, const int ascend);
static void m2startif_d (DBL_D *X, size_t N, const int ascend);
static void m2startif_c (CFLT_F *X, size_t N, const int ascend);
static void m2startif_z (CDBL_D *X, size_t N, const int ascend);


static void m2startif_s (FLT_F *X, size_t N, const int ascend)
{
    size_t i = 0u;
    FLT_F m = X[0];
    ++X;
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X[0].val<m.val) { m = X[0]; i = n; } }
    }
    else
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X[0].val>m.val) { m = X[0]; i = n; } }
    }
    X -= N; X[i] = X[0]; X[0] = m;
}


static void m2startif_d (DBL_D *X, size_t N, const int ascend)
{
    size_t i = 0u;
    DBL_D m = X[0];
    ++X;
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X[0].val<m.val) { m = X[0]; i = n; } }
    }
    else
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X[0].val>m.val) { m = X[0]; i = n; } }
    }
    X -= N; X[i] = X[0]; X[0] = m;
}


static void m2startif_c (CFLT_F *X, size_t N, const int ascend)
{
    size_t i = 0u;
    CFLT_F m = X[0];
    ++X;
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X[0].val<m.val) { m = X[0]; i = n; } }
    }
    else
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X[0].val>m.val) { m = X[0]; i = n; } }
    }
    X -= N; X[i] = X[0]; X[0] = m;
}


static void m2startif_z (CDBL_D *X, size_t N, const int ascend)
{
    size_t i = 0u;
    CDBL_D m = X[0];
    ++X;
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X[0].val<m.val) { m = X[0]; i = n; } }
    }
    else
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X[0].val>m.val) { m = X[0]; i = n; } }
    }
    X -= N; X[i] = X[0]; X[0] = m;
}


void partial_sortif_s (FLT_F *X, size_t N, size_t k, const int ascend)
{
    if (k<K_THRESH)
    {
        for (size_t j=0u; j<=k; ++j, --N, ++X) { m2startif_s(X, N, ascend); }
    }
    else
    {
        if (k<N-1u)
        {
            size_t p, cnt=0u;
            while (N>1u)
            {
                p = lomuto_partitionif_s(X, N-1u, N-1u, !ascend);
                if (k==p) { break; }
                else if (k<p) { N = p; }
                else { ++p; X += p; N -= p; k -= p; cnt += p; }
            }
            X -= cnt; k += cnt;
        }
        quicksortif_s(X, k+1u, ascend);
    }
}


void partial_sortif_d (DBL_D *X, size_t N, size_t k, const int ascend)
{
    if (k<K_THRESH)
    {
        for (size_t j=0u; j<=k; ++j, --N, ++X) { m2startif_d(X, N, ascend); }
    }
    else
    {
        if (k<N-1u)
        {
            size_t p, cnt=0u;
            while (N>1u)
            {
                p = lomuto_partitionif_d(X, N-1u, N-1u, !ascend);
                if (k==p) { break; }
                else if (k<p) { N = p; }
                else { ++p; X += p; N -= p; k -= p; cnt += p; }
            }
            X -= cnt; k += cnt;
        }
        quicksortif_d(X, k+1u, ascend);
    }
}


void partial_sortif_c (CFLT_F *X, size_t N, size_t k, const int ascend)
{
    if (k<K_THRESH)
    {
        for (size_t j=0u; j<=k; ++j, --N, ++X) { m2startif_c(X, N, ascend); }
    }
    else
    {
        if (k<N-1u)
        {
            size_t p, cnt=0u;
            while (N>1u)
            {
                p = lomuto_partitionif_c(X, N-1u, N-1u, !ascend);
                if (k==p) { break; }
                else if (k<p) { N = p; }
                else { ++p; X += p; N -= p; k -= p; cnt += p; }
            }
            X -= cnt; k += cnt;
        }
        quicksortif_c(X, k+1u, ascend);
    }
}


void partial_sortif_z (CDBL_D *X, size_t N, size_t k, const int ascend)
{
    if (k<K_THRESH)
    {
        for (size_t j=0u; j<=k; ++j, --N, ++X) { m2startif_z(X, N, ascend); }
    }
    else
    {
        if (k<N-1u)
        {
            size_t p, cnt=0u;
            while (N>1u)
            {
                p = lomuto_partitionif_z(X, N-1u, N-1u, !ascend);
                if (k==p) { break; }
                else if (k<p) { N = p; }
                else { ++p; X += p; N -= p; k -= p; cnt += p; }
            }
            X -= cnt; k += cnt;
        }
        quicksortif_z(X, k+1u, ascend);
    }
}


#ifdef __cplusplus
}
}
#endif
