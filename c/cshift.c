//Circular shift of X by P elements along dim.
//This has no in-place version (no faster version that I can find).

#include <stdio.h>
#include <stdlib.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int cshift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in cshift_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0u) {}
    else if (D==0)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        X += L2;
        for (size_t l=L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
        X -= L;
        for (size_t l=L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += L2;
            for (size_t v=V; v>0u; --v, X+=L)
            {
                for (size_t l=L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
                X -= L;
                for (size_t l=L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            X += V * L2;
            for (size_t n=V*L1; n>0u; --n, ++X, ++Y) { *Y = *X; }
            X -= N;
            for (size_t n=V*L2; n>0u; --n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                X += K * L2;
                for (size_t b=B; b>0u; --b, X-=K*L2-1u, Y-=K*L-1u)
                {
                    for (size_t l=L1; l>0u; --l, X+=K, Y+=K) { *Y = *X; }
                    X -= K*L;
                    for (size_t l=L2; l>0u; --l, X+=K, Y+=K) { *Y = *X; }
                }
            }
        }
    }
    
    return 0;
}


int cshift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in cshift_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0u) {}
    else if (D==0)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        X += L2;
        for (size_t l=L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
        X -= L;
        for (size_t l=L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += L2;
            for (size_t v=V; v>0u; --v, X+=L)
            {
                for (size_t l=L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
                X -= L;
                for (size_t l=L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            X += V * L2;
            for (size_t n=V*L1; n>0u; --n, ++X, ++Y) { *Y = *X; }
            X -= N;
            for (size_t n=V*L2; n>0u; --n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                X += K * L2;
                for (size_t b=B; b>0u; --b, X-=K*L2-1u, Y-=K*L-1u)
                {
                    for (size_t l=L1; l>0u; --l, X+=K, Y+=K) { *Y = *X; }
                    X -= K*L;
                    for (size_t l=L2; l>0u; --l, X+=K, Y+=K) { *Y = *X; }
                }
            }
        }
    }
    
    return 0;
}


int cshift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in cshift_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0u) {}
    else if (D==0)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        X += 2u*L2;
        for (size_t l=2u*L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
        X -= 2u*L;
        for (size_t l=2u*L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += 2u*L2;
            for (size_t v=V; v>0u; --v, X+=2u*L)
            {
                for (size_t l=2u*L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
                X -= 2u*L;
                for (size_t l=2u*L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            X += 2u*V*L2;
            for (size_t n=2u*V*L1; n>0u; --n, ++X, ++Y) { *Y = *X; }
            X -= 2u*N;
            for (size_t n=2u*V*L2; n>0u; --n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                X += 2u*K*L2;
                for (size_t b=B; b>0u; --b, X-=2u*K*L2-2u, Y-=2u*K*L-2u)
                {
                    for (size_t l=L1; l>0u; --l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
                    X -= 2u*K*L;
                    for (size_t l=L2; l>0u; --l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
                }
            }
        }
    }
    
    return 0;
}


int cshift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in cshift_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0u) {}
    else if (D==0)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        X += 2u*L2;
        for (size_t l=2u*L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
        X -= 2u*L;
        for (size_t l=2u*L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += 2u*L2;
            for (size_t v=V; v>0u; --v, X+=2u*L)
            {
                for (size_t l=2u*L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
                X -= 2u*L;
                for (size_t l=2u*L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            X += 2u*V*L2;
            for (size_t n=2u*V*L1; n>0u; --n, ++X, ++Y) { *Y = *X; }
            X -= 2u*N;
            for (size_t n=2u*V*L2; n>0u; --n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                X += 2u*K*L2;
                for (size_t b=B; b>0u; --b, X-=2u*K*L2-2u, Y-=2u*K*L-2u)
                {
                    for (size_t l=L1; l>0u; --l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
                    X -= 2u*K*L;
                    for (size_t l=L2; l>0u; --l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
                }
            }
        }
    }
    
    return 0;
}


int cshift_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in cshift_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0u || D==0) {}
    else
    {
        float *Y;
        if (!(Y=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in cshift_s: problem with malloc. "); perror("malloc"); return 1; }
        
        if (L==N)
        {
            X += L2;
            for (size_t l=L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
            X -= L;
            for (size_t l=L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
            X -= L2; Y -= L;
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *X = *Y; }
            Y -= L;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y-=L)
                {
                    X += L2;
                    for (size_t l=L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    X -= L;
                    for (size_t l=L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    X -= L2; Y -= L;
                    for (size_t l=L; l>0u; --l, ++X, ++Y) { *X = *Y; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=L)
                    {
                        X += K*L2;
                        for (size_t l=L1; l>0u; --l, X+=K, ++Y) { *Y = *X; }
                        X -= K*L;
                        for (size_t l=L2; l>0u; --l, X+=K, ++Y) { *Y = *X; }
                        X -= K*L2; Y -= L;
                        for (size_t l=L; l>0u; --l, X+=K, ++Y) { *X = *Y; }
                    }
                }
            }
        }
        free(Y);
    }
    
    return 0;
}


int cshift_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in cshift_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0u || D==0) {}
    else
    {
        double *Y;
        if (!(Y=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in cshift_s: problem with malloc. "); perror("malloc"); return 1; }
        
        if (L==N)
        {
            X += L2;
            for (size_t l=L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
            X -= L;
            for (size_t l=L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
            X -= L2; Y -= L;
            for (size_t l=L; l>0u; --l, ++X, ++Y) { *X = *Y; }
            Y -= L;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y-=L)
                {
                    X += L2;
                    for (size_t l=L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    X -= L;
                    for (size_t l=L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    X -= L2; Y -= L;
                    for (size_t l=L; l>0u; --l, ++X, ++Y) { *X = *Y; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=L)
                    {
                        X += K*L2;
                        for (size_t l=L1; l>0u; --l, X+=K, ++Y) { *Y = *X; }
                        X -= K*L;
                        for (size_t l=L2; l>0u; --l, X+=K, ++Y) { *Y = *X; }
                        X -= K*L2; Y -= L;
                        for (size_t l=L; l>0u; --l, X+=K, ++Y) { *X = *Y; }
                    }
                }
            }
        }
        free(Y);
    }
    
    return 0;
}


int cshift_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in cshift_inplace_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0u || D==0) {}
    else
    {
        float *Y;
        if (!(Y=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in cshift_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
        
        if (L==N)
        {
            X += 2u*L2;
            for (size_t l=2u*L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
            X -= 2u*L;
            for (size_t l=2u*L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
            X -= 2u*L2; Y -= 2u*L;
            for (size_t l=2u*L; l>0u; --l, ++X, ++Y) { *X = *Y; }
            Y -= 2u*L;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y-=2u*L)
                {
                    X += 2u*L2;
                    for (size_t l=2u*L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    X -= 2u*L;
                    for (size_t l=2u*L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    X -= 2u*L2; Y -= 2u*L;
                    for (size_t l=2u*L; l>0u; --l, ++X, ++Y) { *X = *Y; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*L2-2u, Y-=2u*L)
                    {
                        X += 2u*K*L2;
                        for (size_t l=L1; l>0u; --l, X+=2u*K-1u, ++Y) { *Y = *X; *++Y = *++X; }
                        X -= 2u*K*L;
                        for (size_t l=L2; l>0u; --l, X+=2u*K-1u, ++Y) { *Y = *X; *++Y = *++X; }
                        X -= 2u*K*L2; Y -= 2u*L;
                        for (size_t l=L; l>0u; --l, X+=2u*K-1u, ++Y) { *X = *Y; *++X = *++Y; }
                    }
                }
            }
        }
        free(Y);
    }
    
    return 0;
}


int cshift_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in cshift_inplace_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0u || D==0) {}
    else
    {
        double *Y;
        if (!(Y=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in cshift_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
        
        if (L==N)
        {
            X += 2u*L2;
            for (size_t l=2u*L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
            X -= 2u*L;
            for (size_t l=2u*L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
            X -= 2u*L2; Y -= 2u*L;
            for (size_t l=2u*L; l>0u; --l, ++X, ++Y) { *X = *Y; }
            Y -= 2u*L;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y-=2u*L)
                {
                    X += 2u*L2;
                    for (size_t l=2u*L1; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    X -= 2u*L;
                    for (size_t l=2u*L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    X -= 2u*L2; Y -= 2u*L;
                    for (size_t l=2u*L; l>0u; --l, ++X, ++Y) { *X = *Y; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*L2-2u, Y-=2u*L)
                    {
                        X += 2u*K*L2;
                        for (size_t l=L1; l>0u; --l, X+=2u*K-1u, ++Y) { *Y = *X; *++Y = *++X; }
                        X -= 2u*K*L;
                        for (size_t l=L2; l>0u; --l, X+=2u*K-1u, ++Y) { *Y = *X; *++Y = *++X; }
                        X -= 2u*K*L2; Y -= 2u*L;
                        for (size_t l=L; l>0u; --l, X+=2u*K-1u, ++Y) { *X = *Y; *++X = *++Y; }
                    }
                }
            }
        }
        free(Y);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
