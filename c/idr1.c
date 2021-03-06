//Vec2vec operation.
//Scales each vector in X along dim such that the IDR (inter-decile range) is from 0 to 1.
//If m1, then scales IDR to -1 to 1.
//This operates in-place.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int idr1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char m1);
int idr1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char m1);


int idr1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char m1)
{
    if (dim>3) { fprintf(stderr,"error in idr1_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (L<2) { fprintf(stderr,"error in idr1_s: L (vec length) must be > 1\n"); return 1; }
    float mn, mx, rng;

    //Prep interpolation
    const float p1 = 0.1f*(L-1), p2 = 0.9f*(L-1);
    const size_t i1 = (size_t)floorf(p1), i2 = (size_t)floorf(p2);
    const float w2 = p1 - floorf(p1), w1 = 1.0f - w2;
    const float w4 = p2 - floorf(p2), w3 = 1.0f - w4;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in idr1_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in idr1_s: problem with LAPACKE function\n"); }
        X1 += i1;
        mn = w1**X1 + w2**(X1+1);
        X1 += i2 - i1;
        mx = w3**X1 + w4**(X1+1);
        X1 -= i2;
        rng = mx - mn;
        if (m1)
        {
            for (size_t l=0u; l<L; ++l) { --X; *X = 2.0f*(*X-mn)/rng - 1.0f; }
        }
        else
        {
            for (size_t l=0u; l<L; ++l) { --X; *X = (*X-mn)/rng; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in idr1_s: problem with LAPACKE function\n"); }
                X1 += i1;
                mn = w1**X1 + w2**(X1+1);
                X1 += i2 - i1;
                mx = w3**X1 + w4**(X1+1);
                X1 -= i2;
                rng = mx - mn;
                if (m1)
                {
                    for (size_t l=0u; l<L; ++l, ++X) { *X = 2.0f*(*X-mn)/rng - 1.0f; }
                }
                else
                {
                    for (size_t l=0u; l<L; ++l, ++X) { *X = (*X-mn)/rng; }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, ++X)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in idr1_s: problem with LAPACKE function\n"); }
                    X1 += i1;
                    mn = w1**X1 + w2**(X1+1);
                    X1 += i2 - i1;
                    mx = w3**X1 + w4**(X1+1);
                    X1 -= i2;
                    rng = mx - mn;
                    if (m1)
                    {
                        for (size_t l=0u; l<L; ++l) { X-=K; *X = 2.0f*(*X-mn)/rng - 1.0f; }
                    }
                    else
                    {
                        for (size_t l=0u; l<L; ++l) { X-=K; *X = (*X-mn)/rng; }
                    }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int idr1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char m1)
{
    if (dim>3) { fprintf(stderr,"error in idr1_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (L<2) { fprintf(stderr,"error in idr1_d: L (vec length) must be > 1\n"); return 1; }
    double mn, mx, rng;

    //Prep interpolation
    const double p1 = 0.1*(L-1), p2 = 0.9*(L-1);
    const size_t i1 = (size_t)floor(p1), i2 = (size_t)floor(p2);
    const double w2 = p1 - floor(p1), w1 = 1.0 - w2;
    const double w4 = p2 - floor(p2), w3 = 1.0 - w4;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in idr1_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in idr1_d: problem with LAPACKE function\n"); }
        X1 += i1;
        mn = w1**X1 + w2**(X1+1);
        X1 += i2 - i1;
        mx = w3**X1 + w4**(X1+1);
        X1 -= i2;
        rng = mx - mn;
        if (m1)
        {
            for (size_t l=0u; l<L; ++l) { --X; *X = 2.0*(*X-mn)/rng - 1.0; }
        }
        else
        {
            for (size_t l=0u; l<L; ++l) { --X; *X = (*X-mn)/rng; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in idr1_d: problem with LAPACKE function\n"); }
                X1 += i1;
                mn = w1**X1 + w2**(X1+1);
                X1 += i2 - i1;
                mx = w3**X1 + w4**(X1+1);
                X1 -= i2;
                rng = mx - mn;
                if (m1)
                {
                    for (size_t l=0u; l<L; ++l, ++X) { *X = 2.0*(*X-mn)/rng - 1.0; }
                }
                else
                {
                    for (size_t l=0u; l<L; ++l, ++X) { *X = (*X-mn)/rng; }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, ++X)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in idr1_d: problem with LAPACKE function\n"); }
                    X1 += i1;
                    mn = w1**X1 + w2**(X1+1);
                    X1 += i2 - i1;
                    mx = w3**X1 + w4**(X1+1);
                    X1 -= i2;
                    rng = mx - mn;
                    if (m1)
                    {
                        for (size_t l=0u; l<L; ++l) { X-=K; *X = 2.0*(*X-mn)/rng - 1.0; }
                    }
                    else
                    {
                        for (size_t l=0u; l<L; ++l) { X-=K; *X = (*X-mn)/rng; }
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
