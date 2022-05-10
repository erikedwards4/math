//Sort_Help function.
//Implements insertion sort algorithm (see insert_sort.c for use).

#pragma once

#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


void insertion_sort_s (float *X, const size_t hi, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            float x = *(X+1);
            while (j>0u && *X>x) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    else
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            float x = *(X+1);
            while (j>0u && *X<x) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    X -= hi - 1u;
}


void insertion_sort_d (double *X, const size_t hi, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            double x = *(X+1);
            while (j>0u && *X>x) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    else
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            double x = *(X+1);
            while (j>0u && *X<x) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    X -= hi - 1u;
}


void insertion_sort_c (float *X, const size_t hi, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            float xr = *(X+2), xi = *(X+3), x = xr*xr + xi*xi + 1e-5f;
            while (j>0u && *X**X+*(X+1)**(X+1)>x)
            {
                *(X+2) = *X; *(X+3) = *(X+1);
                X -= 2; --j;
            }
            X += 2;
            *X = xr; *(X+1) = xi;
            X += 2u*(i-j);
        }
    }
    else
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            float xr = *(X+2), xi = *(X+3), x = xr*xr + xi*xi - 1e-5f;
            while (j>0u && *X**X+*(X+1)**(X+1)<x)
            {
                *(X+2) = *X; *(X+3) = *(X+1);
                X -= 2; --j;
            }
            X += 2;
            *X = xr; *(X+1) = xi;
            X += 2u*(i-j);
        }
    }
    X -= 2u*hi - 2u;
}


void insertion_sort_z (double *X, const size_t hi, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            double xr = *(X+2), xi = *(X+3), x = xr*xr + xi*xi + 1e-9;
            while (j>0u && *X**X+*(X+1)**(X+1)>x)
            {
                *(X+2) = *X; *(X+3) = *(X+1);
                X -= 2; --j;
            }
            X += 2;
            *X = xr; *(X+1) = xi;
            X += 2u*(i-j);
        }
    }
    else
    {
        for (size_t i=1u; i<hi; ++i)
        {
            size_t j = i;
            double xr = *(X+2), xi = *(X+3), x = xr*xr + xi*xi - 1e-9;
            while (j>0u && *X**X+*(X+1)**(X+1)<x)
            {
                *(X+2) = *X; *(X+3) = *(X+1);
                X -= 2; --j;
            }
            X += 2;
            *X = xr; *(X+1) = xi;
            X += 2u*(i-j);
        }
    }
    X -= 2u*hi - 2u;
}


#ifdef __cplusplus
}
}
#endif
