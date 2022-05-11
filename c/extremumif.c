//Sortif/Select Help function.
//Gets extremum of vector X with length N.
//This is the min if min is True, else the max.
//This returns the value and is for internal C use only (like all the Help functions).

#pragma once

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


FLT_F extremumif_s (const FLT_F *X, const size_t N, const int min)
{
    size_t i = 0u;
    float m = X->val;
    ++X;
    if (min)
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X->val<m) { m = X->val; i = n; } }
    }
    else
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X->val>m) { m = X->val; i = n; } }
    }
    X -= N - i;
    return *X;
}


DBL_D extremumif_d (const DBL_D *X, const size_t N, const int min)
{
    size_t i = 0u;
    double m = X->val;
    ++X;
    if (min)
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X->val<m) { m = X->val; i = n; } }
    }
    else
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X->val>m) { m = X->val; i = n; } }
    }
    X -= N - i;
    return *X;
}


CFLT_F extremumif_c (const CFLT_F *X, const size_t N, const int min)
{
    size_t i = 0u;
    float m = X->val;
    ++X;
    if (min)
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X->val<m) { m = X->val; i = n; } }
    }
    else
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X->val>m) { m = X->val; i = n; } }
    }
    X -= N - i;
    return *X;
}


CDBL_D extremumif_z (const CDBL_D *X, const size_t N, const int min)
{
    size_t i = 0u;
    double m = X->val;
    ++X;
    if (min)
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X->val<m) { m = X->val; i = n; } }
    }
    else
    {
        for (size_t n=1u; n<N; ++n, ++X) { if (X->val>m) { m = X->val; i = n; } }
    }
    X -= N - i;
    return *X;
}


#ifdef __cplusplus
}
}
#endif
