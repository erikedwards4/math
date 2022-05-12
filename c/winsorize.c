//Vec2vec operation.
//Winsorizes each vector in X along dim.
//This has in-place and not-in-place versions.

//For each vector, winsorizing works as follows:
//The min value above the pth percentile replaces all values below it.
//The max value below the (1-q)th percentile replaces all values above it.
//The output vector has the same length and overall order as the input vector.

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "codee_math.h"
#include "kselect.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int winsorize_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3u) { fprintf(stderr,"error in winsorize_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>=50.0f) { fprintf(stderr,"error in winsorize_s: p1 must be in [0 50)"); return 1; }
    if (q<0.0f || q>=50.0f) { fprintf(stderr,"error in winsorize_s: p2 must be in [0 50)"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    
    const float p1 = (p/100.0f)*(float)(L-1u), p2 = (1.0f-q/100.0f)*(float)(L-1u);
    const size_t i1 = (size_t)ceilf(p1), i2 = (size_t)floorf(p2);
    float x1, x2;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in winsorize_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u || p<=FLT_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X -= L; X1 -= L;
        x2 = kselect_s(X1,L,i2,1);
        x1 = kselect_s(X1,i2,i1,1);
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = (*X<x1) ? x1 : (*X>x2) ? x2 : *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X -= L; X1 -= L;
                x2 = kselect_s(X1,L,i2,1);
                x1 = kselect_s(X1,i2,i1,1);
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = (*X<x1) ? x1 : (*X>x2) ? x2 : *X; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X -= K*L; X1 -= L;
                    x2 = kselect_s(X1,L,i2,1);
                    x1 = kselect_s(X1,i2,i1,1);
                    for (size_t l=L; l>0u; --l, X+=K, Y+=K) { *Y = (*X<x1) ? x1 : (*X>x2) ? x2 : *X; }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int winsorize_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3u) { fprintf(stderr,"error in winsorize_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in winsorize_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    const double p1 = (p/100.0)*(double)(L-1u), p2 = (1.0-q/100.0)*(double)(L-1u);
    const size_t i1 = (size_t)ceil(p1), i2 = (size_t)floor(p2);
    double x1, x2;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in winsorize_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u || p<=DBL_EPSILON)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X -= L; X1 -= L;
        x2 = kselect_d(X1,L,i2,1);
        x1 = kselect_d(X1,i2,i1,1);
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = (*X<x1) ? x1 : (*X>x2) ? x2 : *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X -= L; X1 -= L;
                x2 = kselect_d(X1,L,i2,1);
                x1 = kselect_d(X1,i2,i1,1);
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = (*X<x1) ? x1 : (*X>x2) ? x2 : *X; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X -= K*L; X1 -= L;
                    x2 = kselect_d(X1,L,i2,1);
                    x1 = kselect_d(X1,i2,i1,1);
                    for (size_t l=L; l>0u; --l, X+=K, Y+=K) { *Y = (*X<x1) ? x1 : (*X>x2) ? x2 : *X; }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int winsorize_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3u) { fprintf(stderr,"error in winsorize_inplace_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in winsorize_inplace_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    const float p1 = (p/100.0f)*(float)(L-1u), p2 = (1.0f-q/100.0f)*(float)(L-1u);
    const size_t i1 = (size_t)ceilf(p1), i2 = (size_t)floorf(p2);
    float x1, x2;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in winsorize_inplace_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u || L==1u || p<=FLT_EPSILON) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X -= L; X1 -= L;
        x2 = kselect_s(X1,L,i2,1);
        x1 = kselect_s(X1,i2,i1,1);
        for (size_t l=L; l>0u; --l, ++X)
        {
            if (*X<x1) { *X = x1; }
            else if (*X>x2) { *X = x2; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X -= L; X1 -= L;
                x2 = kselect_s(X1,L,i2,1);
                x1 = kselect_s(X1,i2,i1,1);
                for (size_t l=L; l>0u; --l, ++X)
                {
                    if (*X<x1) { *X = x1; }
                    else if (*X>x2) { *X = x2; }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X -= K*L; X1 -= L;
                    x2 = kselect_s(X1,L,i2,1);
                    x1 = kselect_s(X1,i2,i1,1);
                    for (size_t l=L; l>0u; --l, X+=K)
                    {
                        if (*X<x1) { *X = x1; }
                        else if (*X>x2) { *X = x2; }
                    }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int winsorize_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3u) { fprintf(stderr,"error in winsorize_inplace_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in winsorize_inplace_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    const double p1 = (p/100.0)*(double)(L-1u), p2 = (1.0-q/100.0)*(double)(L-1u);
    const size_t i1 = (size_t)ceil(p1), i2 = (size_t)floor(p2);
    double x1, x2;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in winsorize_inplace_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u || L==1u || p<=DBL_EPSILON) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X -= L; X1 -= L;
        x2 = kselect_d(X1,L,i2,1);
        x1 = kselect_d(X1,i2,i1,1);
        for (size_t l=L; l>0u; --l, ++X)
        {
            if (*X<x1) { *X = x1; }
            else if (*X>x2) { *X = x2; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X -= L; X1 -= L;
                x2 = kselect_d(X1,L,i2,1);
                x1 = kselect_d(X1,i2,i1,1);
                for (size_t l=L; l>0u; --l, ++X)
                {
                    if (*X<x1) { *X = x1; }
                    else if (*X>x2) { *X = x2; }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X -= K*L; X1 -= L;
                    x2 = kselect_d(X1,L,i2,1);
                    x1 = kselect_d(X1,i2,i1,1);
                    for (size_t l=L; l>0u; --l, X+=K)
                    {
                        if (*X<x1) { *X = x1; }
                        else if (*X>x2) { *X = x2; }
                    }
                }
            }
        }
    }

    free(X1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
