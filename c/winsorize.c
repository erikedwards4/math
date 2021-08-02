//Vec2vec operation.
//Winsorizes each vector in X along dim.
//This has in-place and not-in-place versions.

//For each vector, winsorizing works as follows:
//The min value above the pth percentile replaces all values below it.
//The max value below the (1-q)th percentile replaces all values above it.
//The output vector has the same length and overall order as the input vector.

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int winsorize_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q);
int winsorize_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q);

int winsorize_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q);
int winsorize_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q);


int winsorize_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3u) { fprintf(stderr,"error in winsorize_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>=50.0f) { fprintf(stderr,"error in winsorize_s: p1 must be in [0 50)"); return 1; }
    if (q<0.0f || q>=50.0f) { fprintf(stderr,"error in winsorize_s: p2 must be in [0 50)"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    
    const float p1 = (p/100.0f)*(float)(L-1u), p2 = (1.0f-q/100.0f)*(float)(L-1u);
    const size_t i1 = (size_t)ceilf(p1), i2 = (size_t)floorf(p2);
    float mn, mx;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in winsorize_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u || p<=FLT_EPSILON)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X -= L; X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsorize_s: problem with LAPACKE function\n"); }
        mn = *(X1+i1); mx = *(X1+i2);
        for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = (*X<mn) ? mn : (*X>mx) ? mx : *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X -= L; X1 -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsorize_s: problem with LAPACKE function\n"); }
                mn = *(X1+i1); mx = *(X1+i2);
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = (*X<mn) ? mn : (*X>mx) ? mx : *X; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X -= K*L; X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsorize_s: problem with LAPACKE function\n"); }
                    mn = *(X1+i1); mx = *(X1+i2);
                    for (size_t l=0u; l<L; ++l, X+=K, Y+=K) { *Y = (*X<mn) ? mn : (*X>mx) ? mx : *X; }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int winsorize_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3u) { fprintf(stderr,"error in winsorize_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in winsorize_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    const double p1 = (p/100.0)*(double)(L-1u), p2 = (1.0-q/100.0)*(double)(L-1u);
    const size_t i1 = (size_t)ceil(p1), i2 = (size_t)floor(p2);
    double mn, mx;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in winsorize_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u || p<=DBL_EPSILON)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X -= L; X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsorize_d: problem with LAPACKE function\n"); }
        mn = *(X1+i1); mx = *(X1+i2);
        for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = (*X<mn) ? mn : (*X>mx) ? mx : *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X -= L; X1 -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsorize_d: problem with LAPACKE function\n"); }
                mn = *(X1+i1); mx = *(X1+i2);
                for (size_t l=0u; l<L; ++l, ++X, ++Y) { *Y = (*X<mn) ? mn : (*X>mx) ? mx : *X; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X -= K*L; X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsorize_d: problem with LAPACKE function\n"); }
                    mn = *(X1+i1); mx = *(X1+i2);
                    for (size_t l=0u; l<L; ++l, X+=K, Y+=K) { *Y = (*X<mn) ? mn : (*X>mx) ? mx : *X; }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int winsorize_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p, const float q)
{
    if (dim>3u) { fprintf(stderr,"error in winsorize_inplace_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>50.0f) { fprintf(stderr,"error in winsorize_inplace_s: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    const float p1 = (p/100.0f)*(float)(L-1u), p2 = (1.0f-q/100.0f)*(float)(L-1u);
    const size_t i1 = (size_t)ceilf(p1), i2 = (size_t)floorf(p2);
    float mn, mx;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in winsorize_inplace_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u || L==1u || p<=FLT_EPSILON) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X -= L; X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsorize_inplace_s: problem with LAPACKE function\n"); }
        mn = *(X1+i1); mx = *(X1+i2);
        for (size_t l=0u; l<L; ++l, ++X)
        {
            if (*X<mn) { *X = mn; }
            else if (*X>mx) { *X = mx; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X -= L; X1 -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsorize_inplace_s: problem with LAPACKE function\n"); }
                mn = *(X1+i1); mx = *(X1+i2);
                for (size_t l=0u; l<L; ++l, ++X)
                {
                    if (*X<mn) { *X = mn; }
                    else if (*X>mx) { *X = mx; }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X -= K*L; X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsorize_inplace_s: problem with LAPACKE function\n"); }
                    mn = *(X1+i1); mx = *(X1+i2);
                    for (size_t l=0u; l<L; ++l, X+=K)
                    {
                        if (*X<mn) { *X = mn; }
                        else if (*X>mx) { *X = mx; }
                    }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int winsorize_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p, const double q)
{
    if (dim>3u) { fprintf(stderr,"error in winsorize_inplace_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>50.0) { fprintf(stderr,"error in winsorize_inplace_d: p must be in [0 50]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    const double p1 = (p/100.0)*(double)(L-1u), p2 = (1.0-q/100.0)*(double)(L-1u);
    const size_t i1 = (size_t)ceil(p1), i2 = (size_t)floor(p2);
    double mn, mx;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in winsorize_inplace_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u || L==1u || p<=DBL_EPSILON) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X -= L; X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsorize_inplace_d: problem with LAPACKE function\n"); }
        mn = *(X1+i1); mx = *(X1+i2);
        for (size_t l=0u; l<L; ++l, ++X)
        {
            if (*X<mn) { *X = mn; }
            else if (*X>mx) { *X = mx; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X -= L; X1 -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsorize_inplace_d: problem with LAPACKE function\n"); }
                mn = *(X1+i1); mx = *(X1+i2);
                for (size_t l=0u; l<L; ++l, ++X)
                {
                    if (*X<mn) { *X = mn; }
                    else if (*X>mx) { *X = mx; }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X -= K*L; X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in winsorize_inplace_d: problem with LAPACKE function\n"); }
                    mn = *(X1+i1); mx = *(X1+i2);
                    for (size_t l=0u; l<L; ++l, X+=K)
                    {
                        if (*X<mn) { *X = mn; }
                        else if (*X>mx) { *X = mx; }
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
