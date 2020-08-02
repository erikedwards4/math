//Vec2scalar (reduction) operation.
//Gets trimmed standard deviation for each vector in X along dim.

//The bottom p% and the top q% of values are excluded.
//These are not percentiles, just percentages of data to exclude,
//although they are approximately equal to the percentiles.

//The inplace version still outputs Y, but modifies X during processing.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int trimstd_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q, const char biased);
int trimstd_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q, const char biased);

int trimstd_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q, const char biased);
int trimstd_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q, const char biased);


int trimstd_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in trimstd_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in trimstd_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const float p1 = (p/100.0f)*(L-1), p2 = (1.0f-q/100.0f)*(L-1);
    size_t i1 = (size_t)floorf(p1), i2 = (size_t)ceilf(p2);
    const size_t Lt = i2 - i1 + 1;
    const float den = 1.0f/Lt, den2 = (biased) ? den : 1.0f/(Lt-1);
    float x, sm, sm2;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in trimstd_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        const float z = 0.0f;
        cblas_scopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        cblas_scopy((int)L,X,1,X1,1);
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimstd_s: problem with LAPACKE function\n"); }
        sm = sm2 = 0.0f;
        X1 += i1;
        for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
        *Y = sm * den;
        for (size_t l=i1; l<=i2; ++l) { x = *--X1 - *Y; sm2 += x*x; }
        X1 -= i1;
        *Y = sqrtf(sm2*den2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L, X1-=i1, ++Y)
            {
                cblas_scopy((int)L,X,1,X1,1);
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimstd_s: problem with LAPACKE function\n"); }
                sm = sm2 = 0.0f;
                for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                *Y = sm * den;
                for (size_t l=i1; l<=i2; ++l) { x = *--X1 - *Y; sm2 += x*x; }
                *Y = sqrtf(sm2*den2);
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, X1-=i1, ++Y)
                {
                    cblas_scopy((int)L,X,(int)K,X1,1);
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimstd_s: problem with LAPACKE function\n"); }
                    sm = sm2 = 0.0f;
                    for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm * den;
                    for (size_t l=i1; l<=i2; ++l) { x = *--X1 - *Y; sm2 += x*x; }
                    *Y = sqrtf(sm2*den2);
                }
            }
        }
    }

    free(X1);
    return 0;
}


int trimstd_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in trimstd_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in trimstd_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const double p1 = (p/100.0)*(L-1), p2 = (1.0-q/100.0)*(L-1);
    size_t i1 = (size_t)floor(p1), i2 = (size_t)ceil(p2);
    const size_t Lt = i2 - i1 + 1;
    const double den = 1.0/Lt, den2 = (biased) ? den : 1.0/(Lt-1);
    double x, sm, sm2;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in trimstd_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        const double z = 0.0;
        cblas_dcopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        cblas_dcopy((int)L,X,1,X1,1);
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimstd_d: problem with LAPACKE function\n"); }
        sm = sm2 = 0.0;
        X1 += i1;
        for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
        *Y = sm * den;
        for (size_t l=i1; l<=i2; ++l) { x = *--X1 - *Y; sm2 += x*x; }
        X1 -= i1;
        *Y = sqrt(sm2*den2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L, X1-=i1, ++Y)
            {
                cblas_dcopy((int)L,X,1,X1,1);
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimstd_d: problem with LAPACKE function\n"); }
                sm = sm2 = 0.0;
                for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                *Y = sm * den;
                for (size_t l=i1; l<=i2; ++l) { x = *--X1 - *Y; sm2 += x*x; }
                *Y = sqrt(sm2*den2);
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, X1-=i1, ++Y)
                {
                    cblas_dcopy((int)L,X,(int)K,X1,1);
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimstd_d: problem with LAPACKE function\n"); }
                    sm = sm2 = 0.0;
                    for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm * den;
                    for (size_t l=i1; l<=i2; ++l) { x = *--X1 - *Y; sm2 += x*x; }
                    *Y = sqrt(sm2*den2);
                }
            }
        }
    }

    free(X1);
    return 0;
}


int trimstd_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p, const float q, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in trimstd_inplace_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in trimstd_inplace_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const float p1 = (p/100.0f)*(L-1), p2 = (1.0f-q/100.0f)*(L-1);
    size_t i1 = (size_t)floorf(p1), i2 = (size_t)ceilf(p2);
    const size_t Lt = i2 - i1 + 1;
    const float den = 1.0f/Lt, den2 = (biased) ? den : 1.0f/(Lt-1);
    float x, sm, sm2;

    if (N==0) {}
    else if (L==1)
    {
        const float z = 0.0f;
        cblas_scopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trimstd_inplace_s: problem with LAPACKE function\n"); }
        if (Lt<7000)
        {
            sm = sm2 = 0.0f;
            X += i1;
            for (size_t l=i1; l<=i2; ++l, ++X) { sm += *X; }
            *Y = sm * den;
            for (size_t l=i1; l<=i2; ++l) { x = *--X - *Y; sm2 += x*x; }
            *Y = sqrtf(sm2*den2);
        }
        else
        {
            const float o = 1.0f;
            *Y = cblas_sdot((int)Lt,&X[i1],1,&o,0) * den;
            cblas_saxpy((int)Lt,-*Y,&o,0,&X[i1],1);
            sm2 = cblas_sdot((int)Lt,&X[i1],1,&X[i1],1);
            *Y = sqrtf(sm2*den2);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L-i1, ++Y)
            {
                if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trimstd_inplace_s: problem with LAPACKE function\n"); }
                sm = sm2 = 0.0f;
                for (size_t l=i1; l<=i2; ++l, ++X) { sm += *X; }
                *Y = sm * den;
                for (size_t l=i1; l<=i2; ++l) { x = *--X - *Y; sm2 += x*x; }
                *Y = sqrtf(sm2*den2);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in trimstd_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, X1-=i1, ++Y)
                {
                    cblas_scopy((int)L,X,(int)K,X1,1);
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimstd_inplace_s: problem with LAPACKE function\n"); }
                    sm = sm2 = 0.0f;
                    for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm * den;
                    for (size_t l=i1; l<=i2; ++l) { x = *--X1 - *Y; sm2 += x*x; }
                    *Y = sqrtf(sm2*den2);
                }
            }
            free(X1);
        }
    }

    return 0;
}


int trimstd_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p, const double q, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in trimstd_inplace_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in trimstd_inplace_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    const double p1 = (p/100.0)*(L-1), p2 = (1.0-q/100.0)*(L-1);
    size_t i1 = (size_t)floor(p1), i2 = (size_t)ceil(p2);
    const size_t Lt = i2 - i1 + 1;
    const double den = 1.0/Lt, den2 = (biased) ? den : 1.0/(Lt-1);
    double x, sm, sm2;

    if (N==0) {}
    else if (L==1)
    {
        const double z = 0.0;
        cblas_dcopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trimstd_inplace_d: problem with LAPACKE function\n"); }
        if (Lt<7000)
        {
            sm = sm2 = 0.0;
            X += i1;
            for (size_t l=i1; l<=i2; ++l, ++X) { sm += *X; }
            *Y = sm * den;
            for (size_t l=i1; l<=i2; ++l) { x = *--X - *Y; sm2 += x*x; }
            *Y = sqrt(sm2*den2);
        }
        else
        {
            const double o = 1.0;
            *Y = cblas_ddot((int)Lt,&X[i1],1,&o,0) * den;
            cblas_daxpy((int)Lt,-*Y,&o,0,&X[i1],1);
            sm2 = cblas_ddot((int)Lt,&X[i1],1,&X[i1],1);
            *Y = sqrt(sm2*den2);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L-i1, ++Y)
            {
                if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in trimstd_inplace_d: problem with LAPACKE function\n"); }
                sm = sm2 = 0.0;
                for (size_t l=i1; l<=i2; ++l, ++X) { sm += *X; }
                *Y = sm * den;
                for (size_t l=i1; l<=i2; ++l) { x = *--X - *Y; sm2 += x*x; }
                *Y = sqrt(sm2*den2);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in trimstd_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, X1-=i1, ++Y)
                {
                    cblas_dcopy((int)L,X,(int)K,X1,1);
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in trimstd_inplace_d: problem with LAPACKE function\n"); }
                    sm = sm2 = 0.0;
                    for (size_t l=i1; l<=i2; ++l, ++X1) { sm += *X1; }
                    *Y = sm * den;
                    for (size_t l=i1; l<=i2; ++l) { x = *--X1 - *Y; sm2 += x*x; }
                    *Y = sqrt(sm2*den2);
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
