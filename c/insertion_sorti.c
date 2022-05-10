//Sorti_Help function.
//Implements insertion sort algorithm (see insert_sorti.c for use).

#pragma once

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


void insertion_sorti_s (FLT_I *X, const size_t hi, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            FLT_I x = *(X+1);
            while (j>0u && X->val>x.val) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    else
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            FLT_I x = *(X+1);
            while (j>0u && X->val<x.val) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    X -= hi - 1u;
}


void insertion_sorti_d (DBL_I *X, const size_t hi, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            DBL_I x = *(X+1);
            while (j>0u && X->val>x.val) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    else
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            DBL_I x = *(X+1);
            while (j>0u && X->val<x.val) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    X -= hi - 1u;
}


void insertion_sorti_c (CFLT_I *X, const size_t hi, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            CFLT_I x = *(X+1);
            while (j>0u && X->val>x.val+1e-5f) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    else
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            CFLT_I x = *(X+1);
            while (j>0u && X->val<x.val-1e-5f) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    X -= hi - 1u;
}


void insertion_sorti_z (CDBL_I *X, const size_t hi, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            CDBL_I x = *(X+1);
            while (j>0u && X->val>x.val+1e-9) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    else
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            CDBL_I x = *(X+1);
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
