//Vec2scalar (reduction) operation.
//Gets interquartile range of each vector in X along dim.
//The IQR is the 75th minus the 25th percentile.

//The in-place versions still return the same Y, but modify X during processing.
//However, it turns out to be almost the identical speed for matrices.

#include <stdio.h>
#include <math.h>
// #include <lapacke.h>
#include "codee_math.h"
#include "extremum.c"
#include "kselect.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int iqr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in iqr_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Prep interpolation
    const float p1 = 0.25f*(float)(L-1u), p2 = 0.75f*(float)(L-1u);
    const size_t i1 = (size_t)floorf(p1), i2 = (size_t)floorf(p2);
    const float w2 = p1 - floorf(p1), w1 = 1.0f - w2;
    const float w4 = p2 - floorf(p2), w3 = 1.0f - w4;
    float x1, x2, x3, x4;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in iqr_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        // if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in iqr_s: problem with LAPACKE function\n"); }
        // *Y = w3*X1[i2] + w4*X1[i2+1u] - w1*X1[i1] - w2*X1[i1+1u];
        x4 = kselect_s(X1,L,i2+1u,1);
        x3 = extremum_s(X1,i2+1u,0);
        x2 = kselect_s(X1,i2,i1+1u,1);
        x1 = extremum_s(X1,i1+1u,0);
        *Y = w3*x3 + w4*x4 - w1*x1 - w2*x2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++Y)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                x4 = kselect_s(X1,L,i2+1u,1);
                x3 = extremum_s(X1,i2+1u,0);
                x2 = kselect_s(X1,i2,i1+1u,1);
                x1 = extremum_s(X1,i1+1u,0);
                *Y = w3*x3 + w4*x4 - w1*x1 - w2*x2;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x4 = kselect_s(X1,L,i2+1u,1);
                    x3 = extremum_s(X1,i2+1u,0);
                    x2 = kselect_s(X1,i2,i1+1u,1);
                    x1 = extremum_s(X1,i1+1u,0);
                    *Y = w3*x3 + w4*x4 - w1*x1 - w2*x2;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int iqr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in iqr_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Prep interpolation
    const double p1 = 0.25*(double)(L-1u), p2 = 0.75*(double)(L-1u);
    const size_t i1 = (size_t)floor(p1), i2 = (size_t)floor(p2);
    const double w2 = p1 - floor(p1), w1 = 1.0 - w2;
    const double w4 = p2 - floor(p2), w3 = 1.0 - w4;
    double x1, x2, x3, x4;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in iqr_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        x4 = kselect_d(X1,L,i2+1u,1);
        x3 = extremum_d(X1,i2+1u,0);
        x2 = kselect_d(X1,i2,i1+1u,1);
        x1 = extremum_d(X1,i1+1u,0);
        *Y = w3*x3 + w4*x4 - w1*x1 - w2*x2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++Y)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                x4 = kselect_d(X1,L,i2+1u,1);
                x3 = extremum_d(X1,i2+1u,0);
                x2 = kselect_d(X1,i2,i1+1u,1);
                x1 = extremum_d(X1,i1+1u,0);
                *Y = w3*x3 + w4*x4 - w1*x1 - w2*x2;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x4 = kselect_d(X1,L,i2+1u,1);
                    x3 = extremum_d(X1,i2+1u,0);
                    x2 = kselect_d(X1,i2,i1+1u,1);
                    x1 = extremum_d(X1,i1+1u,0);
                    *Y = w3*x3 + w4*x4 - w1*x1 - w2*x2;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int iqr_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in iqr_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Prep interpolation
    const float p1 = 0.25f*(float)(L-1u), p2 = 0.75f*(float)(L-1u);
    const size_t i1 = (size_t)floorf(p1), i2 = (size_t)floorf(p2);
    const float w2 = p1 - floorf(p1), w1 = 1.0f - w2;
    const float w4 = p2 - floorf(p2), w3 = 1.0f - w4;
    float x1, x2, x3, x4;

    // struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        // if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in iqr_inplace_s: problem with LAPACKE function\n"); }
        // *Y = w3*X[i2] + w4*X[i2+1u] - w1*X[i1] - w2*X[i1+1u];
        x4 = kselect_s(X,L,i2+1u,1);
        x3 = extremum_s(X,i2+1u,0);
        x2 = kselect_s(X,i2,i1+1u,1);
        x1 = extremum_s(X,i1+1u,0);
        *Y = w3*x3 + w4*x4 - w1*x1 - w2*x2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=L, ++Y)
            {
                // if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in iqr_inplace_s: problem with LAPACKE function\n"); }
                // *Y = w3*X[i2] + w4*X[i2+1u] - w1*X[i1] - w2*X[i1+1u];
                x4 = kselect_s(X,L,i2+1u,1);
                x3 = extremum_s(X,i2+1u,0);
                x2 = kselect_s(X,i2,i1+1u,1);
                x1 = extremum_s(X,i1+1u,0);
                *Y = w3*x3 + w4*x4 - w1*x1 - w2*x2;
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in iqr_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    // if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in iqr_inplace_s: problem with LAPACKE function\n"); }
                    // *Y = w3*X[i2] + w4*X[i2+1u] - w1*X[i1] - w2*X[i1+1u];
                    x4 = kselect_s(X1,L,i2+1u,1);
                    x3 = extremum_s(X1,i2+1u,0);
                    x2 = kselect_s(X1,i2,i1+1u,1);
                    x1 = extremum_s(X1,i1+1u,0);
                    *Y = w3*x3 + w4*x4 - w1*x1 - w2*x2;
                }
            }
            free(X1);
        }
    }

    // clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int iqr_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in iqr_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Prep interpolation
    const double p1 = 0.25*(double)(L-1u), p2 = 0.75*(double)(L-1u);
    const size_t i1 = (size_t)floor(p1), i2 = (size_t)floor(p2);
    const double w2 = p1 - floor(p1), w1 = 1.0 - w2;
    const double w4 = p2 - floor(p2), w3 = 1.0 - w4;
    double x1, x2, x3, x4;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        x4 = kselect_d(X,L,i2+1u,1);
        x3 = extremum_d(X,i2+1u,0);
        x2 = kselect_d(X,i2,i1+1u,1);
        x1 = extremum_d(X,i1+1u,0);
        *Y = w3*x3 + w4*x4 - w1*x1 - w2*x2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=L, ++Y)
            {
                x4 = kselect_d(X,L,i2+1u,1);
                x3 = extremum_d(X,i2+1u,0);
                x2 = kselect_d(X,i2,i1+1u,1);
                x1 = extremum_d(X,i1+1u,0);
                *Y = w3*x3 + w4*x4 - w1*x1 - w2*x2;
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in iqr_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x4 = kselect_d(X1,L,i2+1u,1);
                    x3 = extremum_d(X1,i2+1u,0);
                    x2 = kselect_d(X1,i2,i1+1u,1);
                    x1 = extremum_d(X1,i1+1u,0);
                    *Y = w3*x3 + w4*x4 - w1*x1 - w2*x2;
                }
            }
            free(X1);
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
