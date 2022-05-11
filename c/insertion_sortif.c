//Sortif_Help function.
//Implements insertion sort algorithm (see insert_sortif.c for use).

#pragma once

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


void insertion_sortif_s (FLT_F *X, const size_t N, const int ascend)
{
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            FLT_F x = *(X+1);
            while (m>0u && X->val>x.val) { *(X+1) = *X; --X; --m; }
            *++X = x; X += n - m;
        }
    }
    else
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            FLT_F x = *(X+1);
            while (m>0u && X->val<x.val) { *(X+1) = *X; --X; --m; }
            *++X = x; X += n - m;
        }
    }
    X -= N - 1u;
}


void insertion_sortif_d (DBL_D *X, const size_t N, const int ascend)
{
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            DBL_D x = *(X+1);
            while (m>0u && X->val>x.val) { *(X+1) = *X; --X; --m; }
            *++X = x; X += n - m;
        }
    }
    else
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            DBL_D x = *(X+1);
            while (m>0u && X->val<x.val) { *(X+1) = *X; --X; --m; }
            *++X = x; X += n - m;
        }
    }
    X -= N - 1u;
}


void insertion_sortif_c (CFLT_F *X, const size_t N, const int ascend)
{
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            CFLT_F x = *(X+1);
            while (m>0u && X->val>x.val+1e-5f) { *(X+1) = *X; --X; --m; }
            *++X = x; X += n - m;
        }
    }
    else
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            CFLT_F x = *(X+1);
            while (m>0u && X->val<x.val-1e-5f) { *(X+1) = *X; --X; --m; }
            *++X = x; X += n - m;
        }
    }
    X -= N - 1u;
}


void insertion_sortif_z (CDBL_D *X, const size_t N, const int ascend)
{
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            CDBL_D x = *(X+1);
            while (m>0u && X->val>x.val+1e-9) { *(X+1) = *X; --X; --m; }
            *++X = x; X += n - m;
        }
    }
    else
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            CDBL_D x = *(X+1);
            while (m>0u && X->val<x.val-1e-9) { *(X+1) = *X; --X; --m; }
            *++X = x; X += n - m;
        }
    }
    X -= N - 1u;
}


#ifdef __cplusplus
}
}
#endif
