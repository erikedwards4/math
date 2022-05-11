//Sort_Help function.
//Implements insertion sort algorithm (see insert_sort.c for use).

#pragma once

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


void insertion_sort_s (float *X, const size_t N, const int ascend)
{
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            float x = *(X+1);
            while (m>0u && *X>x) { *(X+1) = *X; --X; --m; }
            *++X = x; X += n - m;
        }
    }
    else
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            float x = *(X+1);
            while (m>0u && *X<x) { *(X+1) = *X; --X; --m; }
            *++X = x; X += n - m;
        }
    }
    X -= N - 1u;
}


void insertion_sort_d (double *X, const size_t N, const int ascend)
{
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            double x = *(X+1);
            while (m>0u && *X>x) { *(X+1) = *X; --X; --m; }
            *++X = x; X += n - m;
        }
    }
    else
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            double x = *(X+1);
            while (m>0u && *X<x) { *(X+1) = *X; --X; --m; }
            *++X = x; X += n - m;
        }
    }
    X -= N - 1u;
}


void insertion_sort_c (float *X, const size_t N, const int ascend)
{
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            float xr = *(X+2), xi = *(X+3), x = xr*xr + xi*xi + 1e-5f;
            while (m>0u && *X**X+*(X+1)**(X+1)>x)
            {
                *(X+2) = *X; *(X+3) = *(X+1);
                X -= 2; --m;
            }
            X += 2;
            *X = xr; *(X+1) = xi;
            X += 2u*(n-m);
        }
    }
    else
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            float xr = *(X+2), xi = *(X+3), x = xr*xr + xi*xi - 1e-5f;
            while (m>0u && *X**X+*(X+1)**(X+1)<x)
            {
                *(X+2) = *X; *(X+3) = *(X+1);
                X -= 2; --m;
            }
            X += 2;
            *X = xr; *(X+1) = xi;
            X += 2u*(n-m);
        }
    }
    X -= 2u*N - 2u;
}


void insertion_sort_z (double *X, const size_t N, const int ascend)
{
    if (ascend)
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            double xr = *(X+2), xi = *(X+3), x = xr*xr + xi*xi + 1e-9;
            while (m>0u && *X**X+*(X+1)**(X+1)>x)
            {
                *(X+2) = *X; *(X+3) = *(X+1);
                X -= 2; --m;
            }
            X += 2;
            *X = xr; *(X+1) = xi;
            X += 2u*(n-m);
        }
    }
    else
    {
        for (size_t n=1u; n<N; ++n)
        {
            size_t m = n;
            double xr = *(X+2), xi = *(X+3), x = xr*xr + xi*xi - 1e-9;
            while (m>0u && *X**X+*(X+1)**(X+1)<x)
            {
                *(X+2) = *X; *(X+3) = *(X+1);
                X -= 2; --m;
            }
            X += 2;
            *X = xr; *(X+1) = xi;
            X += 2u*(n-m);
        }
    }
    X -= 2u*N - 2u;
}


#ifdef __cplusplus
}
}
#endif
