//Makes Toeplitz matrix Y from vec X1 (1st col of Y), and from vec X2 (1st row of Y).
//
//This has versions toeplitz1 and toeplitz2 for 1 input (X1)
//or 2 inputs (X1,X2). In toeplitz1, it is implied that X2==X1.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int toeplitz1_s (float *Y, const float *X, const size_t L)
{
    for (size_t l=0u, l2=0u; l<L; ++l, l2=0u, X-=L-2u*l+1u)
    {
        while (l2<l) { *Y = *X; --X; ++Y; ++l2; }
        while (l2<L) { *Y = *X; ++X; ++Y; ++l2; }
    }

    return 0;
}


int toeplitz1_d (double *Y, const double *X, const size_t L)
{
    for (size_t l=0u, l2=0u; l<L; ++l, l2=0u, X-=L-2u*l+1u)
    {
        while (l2<l) { *Y = *X; --X; ++Y; ++l2; }
        while (l2<L) { *Y = *X; ++X; ++Y; ++l2; }
    }

    return 0;
}


int toeplitz1_c (float *Y, const float *X, const size_t L)
{
    for (size_t l=0u, l2=0u; l<L; ++l, l2=0u, X-=2u*L-4u*l+2u)
    {
        while (l2<l) { *Y = *X; *++Y = *(X+1); X-=2; ++Y; ++l2; }
        while (l2<L) { *Y = *X; *++Y = *++X; ++X; ++Y; ++l2; }
    }

    return 0;
}


int toeplitz1_z (double *Y, const double *X, const size_t L)
{
    for (size_t l=0u, l2=0u; l<L; ++l, l2=0u, X-=2u*L-4u*l+2u)
    {
        while (l2<l) { *Y = *X; *++Y = *(X+1); X-=2; ++Y; ++l2; }
        while (l2<L) { *Y = *X; *++Y = *++X; ++X; ++Y; ++l2; }
    }

    return 0;
}


int toeplitz2_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const int iscolmajor)
{
    if (iscolmajor)
    {
        for (size_t c=C; c>0u; --c, X1-=R-c+1u, X2+=c)
        {
            for (size_t l2=0u; l2<c; ++l2, --X2, ++Y) { *Y = *X2; }
            for (size_t l1=c; l1<R; ++l1, ++X1, ++Y) { *Y = *X1; }
        }
    }
    else
    {
        ++X2;
        for (size_t r=R; r>0u; --r, X1+=r+1u, X2-=C-r)
        {
            for (size_t l1=0u; l1<=r; ++l1, --X1, ++Y) { *Y = *X1; }
            for (size_t l2=r+1u; l2<C; ++l2, ++X2, ++Y) { *Y = *X2; }
        }  
    }

    return 0;
}


int toeplitz2_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const int iscolmajor)
{
    if (iscolmajor)
    {
        for (size_t c=C; c>0u; --c, X1-=R-c+1u, X2+=c)
        {
            for (size_t l2=0u; l2<c; ++l2, --X2, ++Y) { *Y = *X2; }
            for (size_t l1=c; l1<R; ++l1, ++X1, ++Y) { *Y = *X1; }
        }
    }
    else
    {
        ++X2;
        for (size_t r=R; r>0u; --r, X1+=r+1u, X2-=C-r)
        {
            for (size_t l1=0u; l1<=r; ++l1, --X1, ++Y) { *Y = *X1; }
            for (size_t l2=r+1u; l2<C; ++l2, ++X2, ++Y) { *Y = *X2; }
        }
    }

    return 0;
}


int toeplitz2_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const int iscolmajor)
{
    if (iscolmajor)
    {
        for (size_t c=C; c>0u; --c, X1-=2u*(R-c+1u), X2+=2u*c)
        {
            for (size_t l2=0u; l2<c; ++l2, X2-=2, ++Y) { *Y = *X2; *++Y = *(X2+1); }
            for (size_t l1=c; l1<R; ++l1, ++X1, ++Y) { *Y = *X1; *++Y = *++X1; }
        }
    }
    else
    {
        X2+=2;
        for (size_t r=R; r>0u; --r, X1+=2u*r+2u, X2-=2u*(C-r))
        {
            for (size_t l1=0u; l1<=r; ++l1, X1-=2, ++Y) { *Y = *X1; *++Y = *(X1+1); }
            for (size_t l2=r+1u; l2<C; ++l2, ++X2, ++Y) { *Y = *X2; *++Y = *++X2; }
        }
    }

    return 0;
}


int toeplitz2_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const int iscolmajor)
{
    if (iscolmajor)
    {
        for (size_t c=C; c>0u; --c, X1-=2u*(R-c+1u), X2+=2u*c)
        {
            for (size_t l2=0u; l2<c; ++l2, X2-=2, ++Y) { *Y = *X2; *++Y = *(X2+1); }
            for (size_t l1=c; l1<R; ++l1, ++X1, ++Y) { *Y = *X1; *++Y = *++X1; }
        }
    }
    else
    {
        X2+=2;
        for (size_t r=R; r>0u; --r, X1+=2u*r+2u, X2-=2u*(C-r))
        {
            for (size_t l1=0u; l1<=r; ++l1, X1-=2, ++Y) { *Y = *X1; *++Y = *(X1+1); }
            for (size_t l2=r+1u; l2<C; ++l2, ++X2, ++Y) { *Y = *X2; *++Y = *++X2; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
