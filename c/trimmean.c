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

int trimmean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q);
int trimmean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q);

int trimmean_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q);
int trimmean_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q);


int trimmean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3) { fprintf(stderr,"error in trimmean_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in trimmean_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const float p1 = (p/100.0f)*(L-1), p2 = (1.0f-q/100.0f)*(L-1);
    size_t i1 = (size_t)floorf(p1), i2 = (size_t)ceilf(p2);
    const float den = 1.0f / (i2-i1+1);
    float sm;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in trimmean_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimmean_s: problem with LAPACKE function\n"); }
        sm = 0.0f; X1 += i1;
        for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
        X1 -= i2 + 1;
        *Y = sm * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X1-=i2+1, ++Y)
            {
                for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimmean_s: problem with LAPACKE function\n"); }
                sm = 0.0f; X1 += i1;
                for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                *Y = sm * den;
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, X1-=i2+1, ++Y)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
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


int trimmean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3) { fprintf(stderr,"error in trimmean_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in trimmean_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const double p1 = (p/100.0)*(L-1), p2 = (1.0-q/100.0)*(L-1);
    size_t i1 = (size_t)floor(p1), i2 = (size_t)ceil(p2);
    const double den = 1.0 / (i2-i1+1);
    double sm;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in trimmean_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimmean_d: problem with LAPACKE function\n"); }
        sm = 0.0; X1 += i1;
        for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
        X1 -= i2 + 1;
        *Y = sm * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X1-=i2+1, ++Y)
            {
                for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimmean_d: problem with LAPACKE function\n"); }
                sm = 0.0; X1 += i1;
                for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                *Y = sm * den;
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, X1-=i2+1, ++Y)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
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


int trimmean_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3) { fprintf(stderr,"error in trimmean_inplace_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in trimmean_inplace_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const float p1 = (p/100.0f)*(L-1), p2 = (1.0f-q/100.0f)*(L-1);
    size_t i1 = (size_t)floorf(p1), i2 = (size_t)ceilf(p2);
    const float den = 1.0f / (i2-i1+1);
    float sm;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
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
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L-i2-1, ++Y)
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
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, X1-=i2+1, ++Y)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
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


int trimmean_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3) { fprintf(stderr,"error in trimmean_inplace_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in trimmean_inplace_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const double p1 = (p/100.0)*(L-1), p2 = (1.0-q/100.0)*(L-1);
    size_t i1 = (size_t)floor(p1), i2 = (size_t)ceil(p2);
    const double den = 1.0 / (i2-i1+1);
    double sm;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
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
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L-i2-1, ++Y)
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
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, X1-=i2+1, ++Y)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
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
