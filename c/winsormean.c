//Vec2scalar (reduction) operation.
//Gets winsorized mean for each vector in X along dim.

//This replaces values outside of the pth and (1-q)th percentiles with
//the min and max remaining values, respectively, and takes mean.

//The inplace version still outputs Y, but modifies X during processing.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <lapacke.h>
#include "codee_math.h"
#include "kselect.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int winsormean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3u) { fprintf(stderr,"error in winsormean_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>=50.0f) { fprintf(stderr,"error in winsormean_s: p must be in [0 50)"); return 1; }
    if (q<0.0f || q>=50.0f) { fprintf(stderr,"error in winsormean_s: q must be in [0 50)"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = (float)L;

    const float p1 = (p/100.0f)*(float)(L-1u), p2 = (1.0f-q/100.0f)*(float)(L-1u);
    const size_t i1 = (size_t)ceilf(p1), i2 = p2>(float)i1 ? (size_t)floorf(p2) : i1;
    float x1, x2, sm;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in winsormean_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        x2 = kselect_s(X1,L,i2,1);
        x1 = kselect_s(X1,i2,i1,1);
        X1 += i1 + 1u;
        sm = (float)(i1+1u)*x1 + (float)(L-i2)*x2;
        for (size_t l=i1+1u; l<i2; ++l, ++X1) { sm += *X1; }
        X1 -= i2-i1+1u;
        *Y = sm / den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X1-=i2-i1+1u, ++Y)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                x2 = kselect_s(X1,L,i2,1);
                x1 = kselect_s(X1,i2,i1,1);
                X1 += i1 + 1u;
                sm = (float)(i1+1u)*x1 + (float)(L-i2)*x2;
                for (size_t l=i1+1u; l<i2; ++l, ++X1) { sm += *X1; }
                *Y = sm / den;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, X1-=i2-i1+1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x2 = kselect_s(X1,L,i2,1);
                    x1 = kselect_s(X1,i2,i1,1);
                    X1 += i1 + 1u;
                    sm = (float)(i1+1u)*x1 + (float)(L-i2)*x2;
                    for (size_t l=i1+1u; l<i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm / den;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int winsormean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3u) { fprintf(stderr,"error in winsormean_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>=50.0) { fprintf(stderr,"error in winsormean_d: p must be in [0 50)"); return 1; }
    if (q<0.0 || q>=50.0) { fprintf(stderr,"error in winsormean_d: q must be in [0 50)"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = (double)L;

    const double p1 = (p/100.0)*(double)(L-1u), p2 = (1.0-q/100.0)*(double)(L-1u);
    const size_t i1 = (size_t)ceil(p1), i2 = p2>(double)i1 ? (size_t)floor(p2) : i1;
    double x1, x2, sm;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in winsormean_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        x2 = kselect_d(X1,L,i2,1);
        x1 = kselect_d(X1,i2,i1,1);
        X1 += i1 + 1u;
        sm = (double)(i1+1u)*x1 + (double)(L-i2)*x2;
        for (size_t l=i1+1u; l<i2; ++l, ++X1) { sm += *X1; }
        X1 -= i2-i1+1u;
        *Y = sm / den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X1-=i2-i1+1u, ++Y)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                x2 = kselect_d(X1,L,i2,1);
                x1 = kselect_d(X1,i2,i1,1);
                X1 += i1 + 1u;
                sm = (double)(i1+1u)*x1 + (double)(L-i2)*x2;
                for (size_t l=i1+1u; l<i2; ++l, ++X1) { sm += *X1; }
                *Y = sm / den;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, X1-=i2-i1+1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x2 = kselect_d(X1,L,i2,1);
                    x1 = kselect_d(X1,i2,i1,1);
                    X1 += i1 + 1u;
                    sm = (double)(i1+1u)*x1 + (double)(L-i2)*x2;
                    for (size_t l=i1+1u; l<i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm / den;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int winsormean_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3u) { fprintf(stderr,"error in winsormean_inplace_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>=50.0f) { fprintf(stderr,"error in winsormean_inplace_s: p must be in [0 50)"); return 1; }
    if (q<0.0f || q>=50.0f) { fprintf(stderr,"error in winsormean_inplace_s: q must be in [0 50)"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = (float)L;

    const float p1 = (p/100.0f)*(float)(L-1u), p2 = (1.0f-q/100.0f)*(float)(L-1u);
    const size_t i1 = (size_t)ceilf(p1), i2 = p2>(float)i1 ? (size_t)floorf(p2) : i1;
    float x1, x2, sm;

    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        //if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in winsormean_inplace_s: problem with LAPACKE function\n"); }
        //partial_sort_s(X,L,i2,1);
        //X += i2; x2 = *X;
        //X -= i2-i1; x1 = *X++;
        x2 = kselect_s(X,L,i2,1);
        x1 = kselect_s(X,i2,i1,1);
        X += i1 + 1u;
        sm = (float)(i1+1u)*x1 + (float)(L-i2)*x2;
        for (size_t l=i1+1u; l<i2; ++l, ++X) { sm += *X; }
        *Y = sm / den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=L-i2, ++Y)
            {
                x2 = kselect_s(X,L,i2,1);
                x1 = kselect_s(X,i2,i1,1);
                X += i1 + 1u;
                sm = (float)(i1+1u)*x1 + (float)(L-i2)*x2;
                for (size_t l=i1+1u; l<i2; ++l, ++X) { sm += *X; }
                *Y = sm / den;
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in winsormean_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, X1-=i2, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x2 = kselect_s(X1,L,i2,1);
                    x1 = kselect_s(X1,i2,i1,1);
                    X1 += i1 + 1u;
                    sm = (float)(i1+1u)*x1 + (float)(L-i2)*x2;
                    for (size_t l=i1+1u; l<i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm / den;
                }
            }
            free(X1);
        }
    }

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int winsormean_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3u) { fprintf(stderr,"error in winsormean_inplace_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>=50.0) { fprintf(stderr,"error in winsormean_inplace_d: p must be in [0 50)"); return 1; }
    if (q<0.0 || q>=50.0) { fprintf(stderr,"error in winsormean_inplace_d: q must be in [0 50)"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = (double)L;

    const double p1 = (p/100.0)*(double)(L-1u), p2 = (1.0-q/100.0)*(double)(L-1u);
    const size_t i1 = (size_t)ceil(p1), i2 = p2>(double)i1 ? (size_t)floor(p2) : i1;
    double x1, x2, sm;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        x2 = kselect_d(X,L,i2,1);
        x1 = kselect_d(X,i2,i1,1);
        X += i1 + 1u;
        sm = (double)(i1+1u)*x1 + (double)(L-i2)*x2;
        for (size_t l=i1+1u; l<i2; ++l, ++X) { sm += *X; }
        *Y = sm / den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=L-i2, ++Y)
            {
                x2 = kselect_d(X,L,i2,1);
                x1 = kselect_d(X,i2,i1,1);
                X += i1 + 1u;
                sm = (double)(i1+1u)*x1 + (double)(L-i2)*x2;
                for (size_t l=i1+1u; l<i2; ++l, ++X) { sm += *X; }
                *Y = sm / den;
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in winsormean_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, X1-=i2, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x2 = kselect_d(X1,L,i2,1);
                    x1 = kselect_d(X1,i2,i1,1);
                    X1 += i1 + 1u;
                    sm = (double)(i1+1u)*x1 + (double)(L-i2)*x2;
                    for (size_t l=i1+1u; l<i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm / den;
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
