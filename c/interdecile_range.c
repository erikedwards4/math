//Vec2scalar (reduction) operation.
//Gets interdecile range for each vector in X along dim.
//The IDR is the 90th minus the 10th percentile.

//The in-place versions still return the same Y, but modify X during processing.
//However, it turns out to be almost the identical speed for matrices.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int interdecile_range_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int interdecile_range_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int interdecile_range_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int interdecile_range_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int interdecile_range_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in interdecile_range_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    const float p1 = 0.1f*(L-1), p2 = 0.9f*(L-1);
    const size_t i1 = (size_t)floorf(p1), i2 = (size_t)floorf(p2);
    const float w2 = floorf(p1) - p1, w1 = -1.0f - w2;
    const float w4 = p2 - floorf(p2), w3 = 1.0f - w4;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in interdecile_range_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0) {}
    else if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        cblas_scopy((int)L,X,1,X1,1);
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in interdecile_range_s: problem with LAPACKE function\n"); }
        X1 += i1;
        *Y = w1**X1 + w2**(X1+1);
        X1 += i2 - i1;
        *Y += w3**X1 + w4**(X1+1);
        X1 -= i2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L, ++Y)
            {
                cblas_scopy((int)L,X,1,X1,1);
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in interdecile_range_s: problem with LAPACKE function\n"); }
                X1 += i1;
                *Y = w1**X1 + w2**(X1+1);
                X1 += i2 - i1;
                *Y += w3**X1 + w4**(X1+1);
                X1 -= i2;
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_scopy((int)L,X,(int)K,X1,1);
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in interdecile_range_s: problem with LAPACKE function\n"); }
                    X1 += i1;
                    *Y = w1**X1 + w2**(X1+1);
                    X1 += i2 - i1;
                    *Y += w3**X1 + w4**(X1+1);
                    X1 -= i2;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int interdecile_range_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in interdecile_range_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    const double p1 = 0.1*(L-1), p2 = 0.9*(L-1);
    const size_t i1 = (size_t)floor(p1), i2 = (size_t)floor(p2);
    const double w2 = floor(p1) - p1, w1 = -1.0 - w2;
    const double w4 = p2 - floor(p2), w3 = 1.0 - w4;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in interdecile_range_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0) {}
    else if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        cblas_dcopy((int)L,X,1,X1,1);
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in interdecile_range_d: problem with LAPACKE function\n"); }
        X1 += i1;
        *Y = w1**X1 + w2**(X1+1);
        X1 += i2 - i1;
        *Y += w3**X1 + w4**(X1+1);
        X1 -= i2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L, ++Y)
            {
                cblas_dcopy((int)L,X,1,X1,1);
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in interdecile_range_d: problem with LAPACKE function\n"); }
                X1 += i1;
                *Y = w1**X1 + w2**(X1+1);
                X1 += i2 - i1;
                *Y += w3**X1 + w4**(X1+1);
                X1 -= i2;
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_dcopy((int)L,X,(int)K,X1,1);
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in interdecile_range_d: problem with LAPACKE function\n"); }
                    X1 += i1;
                    *Y = w1**X1 + w2**(X1+1);
                    X1 += i2 - i1;
                    *Y += w3**X1 + w4**(X1+1);
                    X1 -= i2;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int interdecile_range_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in interdecile_range_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    const float p1 = 0.1f*(L-1), p2 = 0.9f*(L-1);
    const size_t i1 = (size_t)floorf(p1), i2 = (size_t)floorf(p2);
    const float w2 = floorf(p1) - p1, w1 = -1.0f - w2;
    const float w4 = p2 - floorf(p2), w3 = 1.0f - w4;

    if (N==0) {}
    else if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in interdecile_range_inplace_s: problem with LAPACKE function\n"); }
        X += i1;
        *Y = w1**X + w2**(X+1);
        X += i2 - i1;
        *Y += w3**X + w4**(X+1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L-i2, ++Y)
            {
                if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in interdecile_range_inplace_s: problem with LAPACKE function\n"); }
                X += i1;
                *Y = w1**X + w2**(X+1);
                X += i2 - i1;
                *Y += w3**X + w4**(X+1);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in interdecile_range_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_scopy((int)L,X,(int)K,X1,1);
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in interdecile_range_inplace_s: problem with LAPACKE function\n"); }
                    X1 += i1;
                    *Y = w1**X1 + w2**(X1+1);
                    X1 += i2 - i1;
                    *Y += w3**X1 + w4**(X1+1);
                    X1 -= i2;
                }
            }
            free(X1);
        }
    }

    return 0;
}


int interdecile_range_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in interdecile_range_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    const double p1 = 0.1*(L-1), p2 = 0.9*(L-1);
    const size_t i1 = (size_t)floor(p1), i2 = (size_t)floor(p2);
    const double w2 = floor(p1) - p1, w1 = -1.0 - w2;
    const double w4 = p2 - floor(p2), w3 = 1.0 - w4;

    if (N==0) {}
    else if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in interdecile_range_inplace_d: problem with LAPACKE function\n"); }
        X += i1;
        *Y = w1**X + w2**(X+1);
        X += i2 - i1;
        *Y += w3**X + w4**(X+1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L-i2, ++Y)
            {
                if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in interdecile_range_inplace_d: problem with LAPACKE function\n"); }
                X += i1;
                *Y = w1**X + w2**(X+1);
                X += i2 - i1;
                *Y += w3**X + w4**(X+1);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in interdecile_range_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_dcopy((int)L,X,(int)K,X1,1);
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in interdecile_range_inplace_d: problem with LAPACKE function\n"); }
                    X1 += i1;
                    *Y = w1**X1 + w2**(X1+1);
                    X1 += i2 - i1;
                    *Y += w3**X1 + w4**(X1+1);
                    X1 -= i2;
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
