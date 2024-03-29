//Vec2vec operation.
//Scales each vector in X along dim such that the IQR (inter-quartile range) is from 0 to 1.
//If m1, then scales IQR to -1 to 1.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "codee_math.h"
#include "extremum.c"
#include "kselect.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int iqr1_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1)
{
    if (dim>3u) { fprintf(stderr,"error in iqr1_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in iqr1_s: L (vec length) must be > 1\n"); return 1; }
    float x1, x2, x3, x4, mn, mx, rng;

    //Prep interpolation
    const float p1 = 0.25f*(float)(L-1u), p2 = 0.75f*(float)(L-1u);
    const size_t i1 = (size_t)floorf(p1), i2 = (size_t)floorf(p2);
    const float w2 = p1 - floorf(p1), w1 = 1.0f - w2;
    const float w4 = p2 - floorf(p2), w3 = 1.0f - w4;

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        Y -= L; X -= L;
        x4 = kselect_s(Y,L,i2+1u,1);
        x3 = extremum_s(Y,i2+1u,0);
        x2 = kselect_s(Y,i2,i1+1u,1);
        x1 = extremum_s(Y,i1+1u,0);
        mn = w1*x1 + w2*x2;
        mx = w3*x3 + w4*x4;
        rng = mx - mn;
        if (m1)
        {
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = 2.0f*(*X-mn)/rng - 1.0f; }
        }
        else
        {
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = (*X-mn)/rng; }
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
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
                Y -= L; X -= L;
                x4 = kselect_s(Y,L,i2+1u,1);
                x3 = extremum_s(Y,i2+1u,0);
                x2 = kselect_s(Y,i2,i1+1u,1);
                x1 = extremum_s(Y,i1+1u,0);
                mn = w1*x1 + w2*x2;
                mx = w3*x3 + w4*x4;
                rng = mx - mn;
                if (m1)
                {
                    for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = 2.0f*(*X-mn)/rng - 1.0f; }
                }
                else
                {
                    for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = (*X-mn)/rng; }
                }
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in iqr1_s: problem with malloc. "); perror("malloc"); return 1; }
            
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L; X -= K*L;
                    x4 = kselect_s(X1,L,i2+1u,1);
                    x3 = extremum_s(X1,i2+1u,0);
                    x2 = kselect_s(X1,i2,i1+1u,1);
                    x1 = extremum_s(X1,i1+1u,0);
                    mn = w1*x1 + w2*x2;
                    mx = w3*x3 + w4*x4;
                    rng = mx - mn;
                    if (m1)
                    {
                        for (size_t l=L; l>0u; --l, X+=K, Y+=K) { *Y = 2.0f*(*X-mn)/rng - 1.0f; }
                    }
                    else
                    {
                        for (size_t l=L; l>0u; --l, X+=K, Y+=K) { *Y = (*X-mn)/rng; }
                    }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int iqr1_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1)
{
    if (dim>3u) { fprintf(stderr,"error in iqr1_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in iqr1_d: L (vec length) must be > 1\n"); return 1; }
    double x1, x2, x3, x4, mn, mx, rng;

    //Prep interpolation
    const double p1 = 0.25*(double)(L-1u), p2 = 0.75*(double)(L-1u);
    const size_t i1 = (size_t)floor(p1), i2 = (size_t)floor(p2);
    const double w2 = p1 - floor(p1), w1 = 1.0 - w2;
    const double w4 = p2 - floor(p2), w3 = 1.0 - w4;

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        Y -= L; X -= L;
        x4 = kselect_d(Y,L,i2+1u,1);
        x3 = extremum_d(Y,i2+1u,0);
        x2 = kselect_d(Y,i2,i1+1u,1);
        x1 = extremum_d(Y,i1+1u,0);
        mn = w1*x1 + w2*x2;
        mx = w3*x3 + w4*x4;
        rng = mx - mn;
        if (m1)
        {
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = 2.0*(*X-mn)/rng - 1.0; }
        }
        else
        {
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = (*X-mn)/rng; }
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
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
                Y -= L; X -= L;
                x4 = kselect_d(Y,L,i2+1u,1);
                x3 = extremum_d(Y,i2+1u,0);
                x2 = kselect_d(Y,i2,i1+1u,1);
                x1 = extremum_d(Y,i1+1u,0);
                mn = w1*x1 + w2*x2;
                mx = w3*x3 + w4*x4;
                rng = mx - mn;
                if (m1)
                {
                    for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = 2.0*(*X-mn)/rng - 1.0; }
                }
                else
                {
                    for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = (*X-mn)/rng; }
                }
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in iqr1_d: problem with malloc. "); perror("malloc"); return 1; }
            
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L; X -= K*L;
                    x4 = kselect_d(X1,L,i2+1u,1);
                    x3 = extremum_d(X1,i2+1u,0);
                    x2 = kselect_d(X1,i2,i1+1u,1);
                    x1 = extremum_d(X1,i1+1u,0);
                    mn = w1*x1 + w2*x2;
                    mx = w3*x3 + w4*x4;
                    rng = mx - mn;
                    if (m1)
                    {
                        for (size_t l=L; l>0u; --l, X+=K, Y+=K) { *Y = 2.0*(*X-mn)/rng - 1.0; }
                    }
                    else
                    {
                        for (size_t l=L; l>0u; --l, X+=K, Y+=K) { *Y = (*X-mn)/rng; }
                    }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int iqr1_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1)
{
    if (dim>3u) { fprintf(stderr,"error in iqr1_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in iqr1_inplace_s: L (vec length) must be > 1\n"); return 1; }
    float x1, x2, x3, x4, mn, mx, rng;

    //Prep interpolation
    const float p1 = 0.25f*(float)(L-1u), p2 = 0.75f*(float)(L-1u);
    const size_t i1 = (size_t)floorf(p1), i2 = (size_t)floorf(p2);
    const float w2 = p1 - floorf(p1), w1 = 1.0f - w2;
    const float w4 = p2 - floorf(p2), w3 = 1.0f - w4;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in iqr1_inplace_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        x4 = kselect_s(X1,L,i2+1u,1);
        x3 = extremum_s(X1,i2+1u,0);
        x2 = kselect_s(X1,i2,i1+1u,1);
        x1 = extremum_s(X1,i1+1u,0);
        mn = w1*x1 + w2*x2;
        mx = w3*x3 + w4*x4;
        rng = mx - mn;
        if (m1)
        {
            for (size_t l=L; l>0u; --l) { --X; *X = 2.0f*(*X-mn)/rng - 1.0f; }
        }
        else
        {
            for (size_t l=L; l>0u; --l) { --X; *X = (*X-mn)/rng; }
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
                X1 -= L;
                x4 = kselect_s(X1,L,i2+1u,1);
                x3 = extremum_s(X1,i2+1u,0);
                x2 = kselect_s(X1,i2,i1+1u,1);
                x1 = extremum_s(X1,i1+1u,0);
                mn = w1*x1 + w2*x2;
                mx = w3*x3 + w4*x4;
                rng = mx - mn;
                if (m1)
                {
                    for (size_t l=L; l>0u; --l) { --X; *X = 2.0f*(*X-mn)/rng - 1.0f; }
                }
                else
                {
                    for (size_t l=L; l>0u; --l) { --X; *X = (*X-mn)/rng; }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x4 = kselect_s(X1,L,i2+1u,1);
                    x3 = extremum_s(X1,i2+1u,0);
                    x2 = kselect_s(X1,i2,i1+1u,1);
                    x1 = extremum_s(X1,i1+1u,0);
                    mn = w1*x1 + w2*x2;
                    mx = w3*x3 + w4*x4;
                    rng = mx - mn;
                    if (m1)
                    {
                        for (size_t l=L; l>0u; --l) { X-=K; *X = 2.0f*(*X-mn)/rng - 1.0f; }
                    }
                    else
                    {
                        for (size_t l=L; l>0u; --l) { X-=K; *X = (*X-mn)/rng; }
                    }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int iqr1_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1)
{
    if (dim>3u) { fprintf(stderr,"error in iqr1_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in iqr1_inplace_d: L (vec length) must be > 1\n"); return 1; }
    double x1, x2, x3, x4, mn, mx, rng;

    //Prep interpolation
    const double p1 = 0.25*(double)(L-1u), p2 = 0.75*(double)(L-1u);
    const size_t i1 = (size_t)floor(p1), i2 = (size_t)floor(p2);
    const double w2 = p1 - floor(p1), w1 = 1.0 - w2;
    const double w4 = p2 - floor(p2), w3 = 1.0 - w4;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in iqr1_inplace_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        x4 = kselect_d(X1,L,i2+1u,1);
        x3 = extremum_d(X1,i2+1u,0);
        x2 = kselect_d(X1,i2,i1+1u,1);
        x1 = extremum_d(X1,i1+1u,0);
        mn = w1*x1 + w2*x2;
        mx = w3*x3 + w4*x4;
        rng = mx - mn;
        if (m1)
        {
            for (size_t l=L; l>0u; --l) { --X; *X = 2.0*(*X-mn)/rng - 1.0; }
        }
        else
        {
            for (size_t l=L; l>0u; --l) { --X; *X = (*X-mn)/rng; }
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
                X1 -= L;
                x4 = kselect_d(X1,L,i2+1u,1);
                x3 = extremum_d(X1,i2+1u,0);
                x2 = kselect_d(X1,i2,i1+1u,1);
                x1 = extremum_d(X1,i1+1u,0);
                mn = w1*x1 + w2*x2;
                mx = w3*x3 + w4*x4;
                rng = mx - mn;
                if (m1)
                {
                    for (size_t l=L; l>0u; --l) { --X; *X = 2.0*(*X-mn)/rng - 1.0; }
                }
                else
                {
                    for (size_t l=L; l>0u; --l) { --X; *X = (*X-mn)/rng; }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x4 = kselect_d(X1,L,i2+1u,1);
                    x3 = extremum_d(X1,i2+1u,0);
                    x2 = kselect_d(X1,i2,i1+1u,1);
                    x1 = extremum_d(X1,i1+1u,0);
                    mn = w1*x1 + w2*x2;
                    mx = w3*x3 + w4*x4;
                    rng = mx - mn;
                    if (m1)
                    {
                        for (size_t l=L; l>0u; --l) { X-=K; *X = 2.0*(*X-mn)/rng - 1.0; }
                    }
                    else
                    {
                        for (size_t l=L; l>0u; --l) { X-=K; *X = (*X-mn)/rng; }
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
