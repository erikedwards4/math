//Sorti/Select Help function.
//Gets extremum of vector X with length N.
//This is the min if min is True, else the max.
//This returns the value and is for internal C use only (like all the Help functions).

#pragma once

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


FLT_I extremumi_s (const FLT_I *X, const size_t N, const int min)
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


DBL_I extremumi_d (const DBL_I *X, const size_t N, const int min)
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


CFLT_I extremumi_c (const CFLT_I *X, const size_t N, const int min)
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


CDBL_I extremumi_z (const CDBL_I *X, const size_t N, const int min)
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
