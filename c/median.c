//Vec2scalar (reduction) operation.
//Gets median of each vector in X along dim.

//The in-place versions still return the same Y, but modify X during processing.
//However, it turns out to be almost the identical speed for matrices.

#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int median_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int median_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int median_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int median_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int median_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in median_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t i2 = L/2u;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in median_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in median_s: problem with LAPACKE function\n"); }
        X1 += i2;
        *Y = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
        X1 -= i2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in median_s: problem with LAPACKE function\n"); }
                X1 += i2;
                *Y = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
                X1 -= i2;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in median_s: problem with LAPACKE function\n"); }
                    X1 += i2;
                    *Y = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
                    X1 -= i2;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int median_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in median_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t i2 = L/2u;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in median_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in median_d: problem with LAPACKE function\n"); }
        X1 += i2;
        *Y = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
        X1 -= i2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in median_d: problem with LAPACKE function\n"); }
                X1 += i2;
                *Y = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
                X1 -= i2;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in median_d: problem with LAPACKE function\n"); }
                    X1 += i2;
                    *Y = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
                    X1 -= i2;
                }
            }
        }
    }

    free(X1);
    return 0;
}


int median_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in median_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t i2 = L/2u;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in median_inplace_s: problem with LAPACKE function\n"); }
        X += i2;
        *Y = (L%2u) ? *X : 0.5f*(*X + *(X-1));
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, X+=L-i2, ++Y)
            {
                if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in median_inplace_s: problem with LAPACKE function\n"); }
                X += i2;
                *Y = (L%2u) ? *X : 0.5f*(*X + *(X-1));
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in median_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in median_inplace_s: problem with LAPACKE function\n"); }
                    X1 += i2;
                    *Y = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
                    X1 -= i2;
                }
            }
            free(X1);
        }
    }

    return 0;
}


int median_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in median_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t i2 = L/2u;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in median_inplace_d: problem with LAPACKE function\n"); }
        X += i2;
        *Y = (L%2u) ? *X : 0.5*(*X + *(X-1));
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, X+=L-i2, ++Y)
            {
                if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in median_inplace_d: problem with LAPACKE function\n"); }
                X += i2;
                *Y = (L%2u) ? *X : 0.5*(*X + *(X-1));
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in median_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in median_inplace_d: problem with LAPACKE function\n"); }
                    X1 += i2;
                    *Y = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
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
