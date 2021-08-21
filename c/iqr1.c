//Vec2vec operation.
//Scales each vector in X along dim such that the IQR (inter-quartile range) is from 0 to 1.
//If m1, then scales IQR to -1 to 1.
//This operates in-place.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int iqr1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1);
int iqr1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1);


int iqr1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1)
{
    if (dim>3u) { fprintf(stderr,"error in iqr1_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in iqr1_s: L (vec length) must be > 1\n"); return 1; }
    float mn, mx, rng;

    //Prep interpolation
    const float p1 = 0.25f*(float)(L-1u), p2 = 0.75f*(float)(L-1u);
    const size_t i1 = (size_t)floorf(p1), i2 = (size_t)floorf(p2);
    const float w2 = p1 - floorf(p1), w1 = 1.0f - w2;
    const float w4 = p2 - floorf(p2), w3 = 1.0f - w4;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in iqr1_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in iqr1_s: problem with LAPACKE function\n"); }
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
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in iqr1_s: problem with LAPACKE function\n"); }
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
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in iqr1_s: problem with LAPACKE function\n"); }
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


int iqr1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int m1)
{
    if (dim>3u) { fprintf(stderr,"error in iqr1_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in iqr1_d: L (vec length) must be > 1\n"); return 1; }
    double mn, mx, rng;

    //Prep interpolation
    const double p1 = 0.25*(double)(L-1u), p2 = 0.75*(double)(L-1u);
    const size_t i1 = (size_t)floor(p1), i2 = (size_t)floor(p2);
    const double w2 = p1 - floor(p1), w1 = 1.0 - w2;
    const double w4 = p2 - floor(p2), w3 = 1.0 - w4;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in iqr1_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in iqr1_d: problem with LAPACKE function\n"); }
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
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in iqr1_d: problem with LAPACKE function\n"); }
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
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in iqr1_d: problem with LAPACKE function\n"); }
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
