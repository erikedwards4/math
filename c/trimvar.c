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

int trimvar_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q, const char biased);
int trimvar_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q, const char biased);

int trimvar_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q, const char biased);
int trimvar_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q, const char biased);


int trimvar_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in trimvar_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in trimvar_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const float p1 = (p/100.0f)*(L-1), p2 = (1.0f-q/100.0f)*(L-1);
    size_t i1 = (size_t)floorf(p1), i2 = (size_t)ceilf(p2);
    const size_t Lt = i2 - i1 + 1;
    const float den = 1.0f/Lt, den2 = (biased) ? den : 1.0f/(Lt-1);
    float x, mn, sm2;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in trimvar_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimvar_s: problem with LAPACKE function\n"); }
        mn = sm2 = 0.0f;
        X1 += i1;
        for (size_t l=i1; l<=i2; ++l, ++X1) { mn += *X1; }
        mn *= den;
        for (size_t l=i1; l<=i2; ++l) { x = *--X1 - mn; sm2 += x*x; }
        X1 -= i1;
        *Y = sm2 * den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, X1-=i1, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimvar_s: problem with LAPACKE function\n"); }
                mn = sm2 = 0.0f;
                X1 += i1;
                for (size_t l=i1; l<=i2; ++l, ++X1) { mn += *X1; }
                mn *= den;
                for (size_t l=i1; l<=i2; ++l) { x = *--X1 - mn; sm2 += x*x; }
                *Y = sm2 * den2;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, X1-=i1, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimvar_s: problem with LAPACKE function\n"); }
                    mn = sm2 = 0.0f;
                    X1 += i1;
                    for (size_t l=i1; l<=i2; ++l, ++X1) { mn += *X1; }
                    mn *= den;
                    for (size_t l=i1; l<=i2; ++l) { x = *--X1 - mn; sm2 += x*x; }
                    *Y = sm2 * den2;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int trimvar_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in trimvar_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in trimvar_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const double p1 = (p/100.0)*(L-1), p2 = (1.0-q/100.0)*(L-1);
    size_t i1 = (size_t)floor(p1), i2 = (size_t)ceil(p2);
    const size_t Lt = i2 - i1 + 1;
    const double den = 1.0/Lt, den2 = (biased) ? den : 1.0/(Lt-1);
    double x, mn, sm2;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in trimvar_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimvar_d: problem with LAPACKE function\n"); }
        mn = sm2 = 0.0;
        X1 += i1;
        for (size_t l=i1; l<=i2; ++l, ++X1) { mn += *X1; }
        mn *= den;
        for (size_t l=i1; l<=i2; ++l) { x = *--X1 - mn; sm2 += x*x; }
        X1 -= i1;
        *Y = sm2 * den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, X1-=i1, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimvar_d: problem with LAPACKE function\n"); }
                mn = sm2 = 0.0;
                X1 += i1;
                for (size_t l=i1; l<=i2; ++l, ++X1) { mn += *X1; }
                mn *= den;
                for (size_t l=i1; l<=i2; ++l) { x = *--X1 - mn; sm2 += x*x; }
                *Y = sm2 * den2;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, X1-=i1, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimvar_d: problem with LAPACKE function\n"); }
                    mn = sm2 = 0.0;
                    X1 += i1;
                    for (size_t l=i1; l<=i2; ++l, ++X1) { mn += *X1; }
                    mn *= den;
                    for (size_t l=i1; l<=i2; ++l) { x = *--X1 - mn; sm2 += x*x; }
                    *Y = sm2 * den2;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int trimvar_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in trimvar_inplace_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in trimvar_inplace_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const float p1 = (p/100.0f)*(L-1), p2 = (1.0f-q/100.0f)*(L-1);
    size_t i1 = (size_t)floorf(p1), i2 = (size_t)ceilf(p2);
    const size_t Lt = i2 - i1 + 1;
    const float den = 1.0f/Lt, den2 = (biased) ? den : 1.0f/(Lt-1);
    float x, mn, sm2;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trimvar_inplace_s: problem with LAPACKE function\n"); }
        mn = sm2 = 0.0f;
        X += i1;
        for (size_t l=i1; l<=i2; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=i1; l<=i2; ++l) { x = *--X - mn; sm2 += x*x; }
        *Y = sm2 * den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, X+=L-i1, ++Y)
            {
                if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trimvar_inplace_s: problem with LAPACKE function\n"); }
                mn = sm2 = 0.0f;
                X += i1;
                for (size_t l=i1; l<=i2; ++l, ++X) { mn += *X; }
                mn *= den;
                for (size_t l=i1; l<=i2; ++l) { x = *--X - mn; sm2 += x*x; }
                *Y = sm2 * den2;
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in trimvar_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, X1-=i1, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimvar_inplace_s: problem with LAPACKE function\n"); }
                    mn = sm2 = 0.0f;
                    X1 += i1;
                    for (size_t l=i1; l<=i2; ++l, ++X1) { mn += *X1; }
                    mn *= den;
                    for (size_t l=i1; l<=i2; ++l) { x = *--X1 - mn; sm2 += x*x; }
                    *Y = sm2 * den2;
                }
            }
            free(X1);
        }
    }

    return 0;
}


int trimvar_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in trimvar_inplace_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in trimvar_inplace_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const double p1 = (p/100.0)*(L-1), p2 = (1.0-q/100.0)*(L-1);
    size_t i1 = (size_t)floor(p1), i2 = (size_t)ceil(p2);
    const size_t Lt = i2 - i1 + 1;
    const double den = 1.0/Lt, den2 = (biased) ? den : 1.0/(Lt-1);
    double x, mn, sm2;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trimvar_inplace_d: problem with LAPACKE function\n"); }
        mn = sm2 = 0.0;
        X += i1;
        for (size_t l=i1; l<=i2; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=i1; l<=i2; ++l) { x = *--X - mn; sm2 += x*x; }
        *Y = sm2 * den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, X+=L-i1, ++Y)
            {
                if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trimvar_inplace_d: problem with LAPACKE function\n"); }
                mn = sm2 = 0.0;
                X += i1;
                for (size_t l=i1; l<=i2; ++l, ++X) { mn += *X; }
                mn *= den;
                for (size_t l=i1; l<=i2; ++l) { x = *--X - mn; sm2 += x*x; }
                *Y = sm2 * den2;
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in trimvar_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, X1-=i1, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimvar_inplace_d: problem with LAPACKE function\n"); }
                    mn = sm2 = 0.0;
                    X1 += i1;
                    for (size_t l=i1; l<=i2; ++l, ++X1) { mn += *X1; }
                    mn *= den;
                    for (size_t l=i1; l<=i2; ++l) { x = *--X1 - mn; sm2 += x*x; }
                    *Y = sm2 * den2;
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
