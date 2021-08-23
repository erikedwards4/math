//Vec2scalar (reduction) operation.
//Gets trimmed mean for each vector in X along dim.

//The bottom p% and the top q% of values are excluded.
//These are not percentiles, just percentages of data to exclude,
//although they are approximately equal to the percentiles.

//The inplace version still outputs Y, but modifies X during processing.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int trimmean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q);
int trimmean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q);

int trimmean_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q);
int trimmean_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q);


int trimmean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3u) { fprintf(stderr,"error in trimmean_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in trimmean_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    const float p1 = (p/100.0f)*(float)(L-1u), p2 = (1.0f-q/100.0f)*(float)(L-1u);
    size_t i1 = (size_t)floorf(p1), i2 = (size_t)ceilf(p2);
    const float den = 1.0f / (float)(i2-i1+1u);
    float sm;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in trimmean_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimmean_s: problem with LAPACKE function\n"); }
        sm = 0.0f; X1 += i1;
        for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
        X1 -= i2 + 1u;
        *Y = sm * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X1-=i2+1u, ++Y)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimmean_s: problem with LAPACKE function\n"); }
                sm = 0.0f; X1 += i1;
                for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                *Y = sm * den;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, X1-=i2+1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimmean_s: problem with LAPACKE function\n"); }
                    sm = 0.0f; X1 += i1;
                    for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm * den;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int trimmean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3u) { fprintf(stderr,"error in trimmean_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in trimmean_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    const double p1 = (p/100.0)*(double)(L-1u), p2 = (1.0-q/100.0)*(double)(L-1u);
    size_t i1 = (size_t)floor(p1), i2 = (size_t)ceil(p2);
    const double den = 1.0 / (double)(i2-i1+1u);
    double sm;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in trimmean_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimmean_d: problem with LAPACKE function\n"); }
        sm = 0.0; X1 += i1;
        for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
        X1 -= i2 + 1u;
        *Y = sm * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X1-=i2+1u, ++Y)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimmean_d: problem with LAPACKE function\n"); }
                sm = 0.0; X1 += i1;
                for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                *Y = sm * den;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, X1-=i2+1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimmean_d: problem with LAPACKE function\n"); }
                    sm = 0.0; X1 += i1;
                    for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm * den;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int trimmean_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3u) { fprintf(stderr,"error in trimmean_inplace_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in trimmean_inplace_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    const float p1 = (p/100.0f)*(float)(L-1u), p2 = (1.0f-q/100.0f)*(float)(L-1u);
    size_t i1 = (size_t)floorf(p1), i2 = (size_t)ceilf(p2);
    const float den = 1.0f / (float)(i2-i1+1u);
    float sm;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trimmean_s: problem with LAPACKE function\n"); }
        sm = 0.0f; X += i1;
        for (size_t l=i1; l<=i2; ++l, ++X) { sm += *X; }
        *Y = sm * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=L-i2-1u, ++Y)
            {
                if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trimmean_s: problem with LAPACKE function\n"); }
                sm = 0.0f; X += i1;
                for (size_t l=i1; l<=i2; ++l, ++X) { sm += *X; }
                *Y = sm * den;
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in trimmean_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, X1-=i2+1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimmean_s: problem with LAPACKE function\n"); }
                    sm = 0.0f; X1 += i1;
                    for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm * den;
                }
            }
            free(X1);
        }
    }

    return 0;
}


int trimmean_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3u) { fprintf(stderr,"error in trimmean_inplace_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in trimmean_inplace_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    const double p1 = (p/100.0)*(double)(L-1u), p2 = (1.0-q/100.0)*(double)(L-1u);
    size_t i1 = (size_t)floor(p1), i2 = (size_t)ceil(p2);
    const double den = 1.0 / (double)(i2-i1+1u);
    double sm;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trimmean_d: problem with LAPACKE function\n"); }
        sm = 0.0; X += i1;
        for (size_t l=i1; l<=i2; ++l, ++X) { sm += *X; }
        *Y = sm * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=L-i2-1u, ++Y)
            {
                if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trimmean_d: problem with LAPACKE function\n"); }
                sm = 0.0; X += i1;
                for (size_t l=i1; l<=i2; ++l, ++X) { sm += *X; }
                *Y = sm * den;
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in trimmean_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, X1-=i2+1u, ++Y)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimmean_d: problem with LAPACKE function\n"); }
                    sm = 0.0; X1 += i1;
                    for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm * den;
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
