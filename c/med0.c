//Zeros (subtracts) the median of each vector in X along dim.
//For each vector, estimates the median and then subtracts it from each element.
//This operates in-place.

#include <stdio.h>
#include <stdlib.h>
#include "codee_math.h"
#include "extremum.c"
#include "kselect.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int med0_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in med0_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t i2 = L/2u;
    float med, x2;
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        Y -= L;
        x2 = kselect_s(Y,L,i2,1);
        med = (L%2u) ? x2 : 0.5f*(x2+extremum_s(Y,i2,0));
        for (size_t l=L; l>0u; --l, ++Y) { *Y -= med; }
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
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
                Y -= L;
                x2 = kselect_s(Y,L,i2,1);
                med = (L%2u) ? x2 : 0.5f*(x2+extremum_s(Y,i2,0));
                for (size_t l=L; l>0u; --l, ++Y) { *Y -= med; }
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in med0_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L; X -= K*L;
                    x2 = kselect_s(X1,L,i2,1);
                    med = (L%2u) ? x2 : 0.5f*(x2+extremum_s(X1,i2,0));
                    for (size_t l=L; l>0u; --l, X+=K, Y+=K) { *Y = *X - med; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int med0_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in med0_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t i2 = L/2u;
    double med, x2;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        Y -= L;
        x2 = kselect_d(Y,L,i2,1);
        med = (L%2u) ? x2 : 0.5*(x2+extremum_d(Y,i2,0));
        for (size_t l=L; l>0u; --l, ++Y) { *Y -= med; }
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
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
                Y -= L;
                x2 = kselect_d(Y,L,i2,1);
                med = (L%2u) ? x2 : 0.5*(x2+extremum_d(Y,i2,0));
                for (size_t l=L; l>0u; --l, ++Y) { *Y -= med; }
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in med0_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L; X -= K*L;
                    x2 = kselect_d(X1,L,i2,1);
                    med = (L%2u) ? x2 : 0.5*(x2+extremum_d(X1,i2,0));
                    for (size_t l=L; l>0u; --l, X+=K, Y+=K) { *Y = *X - med; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int med0_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in med0_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t i2 = L/2u;
    float med, x2;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in med0_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        x2 = kselect_s(X1,L,i2,1);
        med = (L%2u) ? x2 : 0.5f*(x2+extremum_s(X1,i2,0));
        for (size_t l=L; l>0u; --l) { --X; *X -= med; }
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
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                x2 = kselect_s(X1,L,i2,1);
                med = (L%2u) ? x2 : 0.5f*(x2+extremum_s(X1,i2,0));
                for (size_t l=L; l>0u; --l, ++X) { *X -= med; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x2 = kselect_s(X1,L,i2,1);
                    med = (L%2u) ? x2 : 0.5f*(x2+extremum_s(X1,i2,0));
                    for (size_t l=L; l>0u; --l) { X-=K; *X -= med; }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int med0_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in med0_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t i2 = L/2u;
    double med, x2;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in med0_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L;
        x2 = kselect_d(X1,L,i2,1);
        med = (L%2u) ? x2 : 0.5*(x2+extremum_d(X1,i2,0));
        for (size_t l=L; l>0u; --l) { --X; *X -= med; }
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
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                x2 = kselect_d(X1,L,i2,1);
                med = (L%2u) ? x2 : 0.5*(x2+extremum_d(X1,i2,0));
                for (size_t l=L; l>0u; --l, ++X) { *X -= med; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    x2 = kselect_d(X1,L,i2,1);
                    med = (L%2u) ? x2 : 0.5*(x2+extremum_d(X1,i2,0));
                    for (size_t l=L; l>0u; --l) { X-=K; *X -= med; }
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
