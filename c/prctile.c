//Vec2scalar (reduction) operation.
//Gets pth percentile for each vector in X along dim.

//The in-place versions still return the same Y, but modify X during processing.
//However, it turns out to be almost the identical speed for matrices.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int prctile_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p);
int prctile_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p);

int prctile_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p);
int prctile_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p);


int prctile_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p)
{
    if (dim>3) { fprintf(stderr,"error in prctile_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>100.0f) { fprintf(stderr,"error in prctile_s: p must be in [0 100]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    const float p1 = (p/100.0f)*(L-1);
    const size_t i1 = (p<100.0f) ? (size_t)floorf(p1) : L-2;
    const float w2 = (p<100.0f) ? p1-floorf(p1) : 1.0f;
    const float w1 = 1.0f - w2;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in prctile_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in prctile_s: problem with LAPACKE function\n"); }
        X1 += i1;
        *Y = w1**X1 + w2**(X1+1);
        X1 -= i1;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in prctile_s: problem with LAPACKE function\n"); }
                X1 += i1;
                *Y = w1**X1 + w2**(X1+1);
                X1 -= i1;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in prctile_s: problem with LAPACKE function\n"); }
                    X1 += i1;
                    *Y = w1**X1 + w2**(X1+1);
                    X1 -= i1;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int prctile_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p)
{
    if (dim>3) { fprintf(stderr,"error in prctile_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>100.0) { fprintf(stderr,"error in prctile_d: p must be in [0 100]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    const double p1 = (p/100.0)*(L-1);
    const size_t i1 = (p<100.0) ? (size_t)floor(p1) : L-2;
    const double w2 = (p<100.0) ? p1-floor(p1) : 1.0;
    const double w1 = 1.0 - w2;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in prctile_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in prctile_d: problem with LAPACKE function\n"); }
        X1 += i1;
        *Y = w1**X1 + w2**(X1+1);
        X1 -= i1;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in prctile_d: problem with LAPACKE function\n"); }
                X1 += i1;
                *Y = w1**X1 + w2**(X1+1);
                X1 -= i1;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in prctile_d: problem with LAPACKE function\n"); }
                    X1 += i1;
                    *Y = w1**X1 + w2**(X1+1);
                    X1 -= i1;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int prctile_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p)
{
    if (dim>3) { fprintf(stderr,"error in prctile_inplace_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>100.0f) { fprintf(stderr,"error in prctile_inplace_s: p must be in [0 100]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    const float p1 = (p/100.0f)*(L-1);
    const size_t i1 = (p<100.0f) ? (size_t)floorf(p1) : L-2;
    const float w2 = (p<100.0f) ? p1-floorf(p1) : 1.0f;
    const float w1 = 1.0f - w2;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in prctile_inplace_s: problem with LAPACKE function\n"); }
        X += i1;
        *Y = w1**X + w2**(X+1);
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
                if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in prctile_inplace_s: problem with LAPACKE function\n"); }
                X += i1;
                *Y = w1**X + w2**(X+1);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in prctile_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in prctile_inplace_s: problem with LAPACKE function\n"); }
                    X1 += i1;
                    *Y = w1**X1 + w2**(X1+1);
                    X1 -= i1;
                }
            }
            free(X1);
        }
    }

    return 0;
}


int prctile_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p)
{
    if (dim>3) { fprintf(stderr,"error in prctile_inplace_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>100.0) { fprintf(stderr,"error in prctile_inplace_d: p must be in [0 100]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    const double p1 = (p/100.0)*(L-1);
    const size_t i1 = (p<100.0) ? (size_t)floor(p1) : L-2;
    const double w2 = (p<100.0) ? p1-floor(p1) : 1.0;
    const double w1 = 1.0 - w2;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in prctile_inplace_d: problem with LAPACKE function\n"); }
        X += i1;
        *Y = w1**X + w2**(X+1);
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
                if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in prctile_inplace_d: problem with LAPACKE function\n"); }
                X += i1;
                *Y = w1**X + w2**(X+1);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in prctile_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in prctile_inplace_d: problem with LAPACKE function\n"); }
                    X1 += i1;
                    *Y = w1**X1 + w2**(X1+1);
                    X1 -= i1;
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
