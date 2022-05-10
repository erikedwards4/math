//Sortif_Help function.
//Implements insertion sort algorithm (see insert_sortif.c for use).

#pragma once

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


void insertion_sortif_s (FLT_F *X, const size_t hi, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            FLT_F x = *(X+1);
            while (j>0u && X->val>x.val) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    else
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            FLT_F x = *(X+1);
            while (j>0u && X->val<x.val) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    X -= hi - 1u;
}


void insertion_sortif_d (DBL_D *X, const size_t hi, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            DBL_D x = *(X+1);
            while (j>0u && X->val>x.val) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    else
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            DBL_D x = *(X+1);
            while (j>0u && X->val<x.val) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    X -= hi - 1u;
}


void insertion_sortif_c (CFLT_F *X, const size_t hi, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            CFLT_F x = *(X+1);
            while (j>0u && X->val>x.val+1e-5f) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    else
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            CFLT_F x = *(X+1);
            while (j>0u && X->val<x.val-1e-5f) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    X -= hi - 1u;
}


void insertion_sortif_z (CDBL_D *X, const size_t hi, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            CDBL_D x = *(X+1);
            while (j>0u && X->val>x.val+1e-9) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    else
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            CDBL_D x = *(X+1);
            while (j>0u && X->val<x.val-1e-9) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    X -= hi - 1u;
}


#ifdef __cplusplus
}
}
#endif
