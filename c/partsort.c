//Vec2vec operation.
//Partial sorts each vector in X along dim using partial_sort (based on quickselect algo).
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <stdlib.h>
#include "codee_math.h"
#include "partial_sort.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int partsort_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend, const size_t k)
{
    if (dim>3u) { fprintf(stderr,"error in partsort_s: dim must be in [0 3]\n"); return 1; }
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (k>L-1u) { fprintf(stderr,"error in partsort_s: k must be in [0 L-1]\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        Y -= L;
        partial_sort_s(Y,L,k,ascend);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, Y+=L)
            {
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
                Y -= L;
                partial_sort_s(Y,L,k,ascend);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in partsort_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, X1-=L, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    partial_sort_s(X1,L,k,ascend);
                    for (size_t l=L; l>0u; --l, ++X1, Y+=K) { *Y = *X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int partsort_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend, const size_t k)
{
    if (dim>3u) { fprintf(stderr,"error in partsort_d: dim must be in [0 3]\n"); return 1; }
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (k>L-1u) { fprintf(stderr,"error in partsort_d: k must be in [0 L-1]\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        partial_sort_d(Y,L,k,ascend);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, Y+=L)
            {
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
                Y -= L;
                partial_sort_d(Y,L,k,ascend);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in partsort_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, X1-=L, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    partial_sort_d(X1,L,k,ascend);
                    for (size_t l=L; l>0u; --l, ++X1, Y+=K) { *Y = *X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int partsort_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend, const size_t k)
{
    if (dim>3u) { fprintf(stderr,"error in partsort_c: dim must be in [0 3]\n"); return 1; }
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (k>L-1u) { fprintf(stderr,"error in partsort_c: k must be in [0 L-1]\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=2u*L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        partial_sort_c(Y,L,k,ascend);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, Y+=2u*L)
            {
                for (size_t l=2u*L; l>0u; --l, ++X, ++Y) { *Y = *X; }
                Y -= 2u*L;
                partial_sort_c(Y,L,k,ascend);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in partsort_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X1-=2u*L, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X1 -= 2u*L;
                    partial_sort_c(X1,L,k,ascend);
                    for (size_t l=L; l>0u; --l, ++X1, Y+=2u*K-1u) { *Y = *X1; *++Y = *++X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int partsort_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend, const size_t k)
{
    if (dim>3u) { fprintf(stderr,"error in partsort_z: dim must be in [0 3]\n"); return 1; }
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (k>L-1u) { fprintf(stderr,"error in partsort_z: k must be in [0 L-1]\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=2u*L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        partial_sort_z(Y,L,k,ascend);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, Y+=2u*L)
            {
                for (size_t l=2u*L; l>0u; --l, ++X, ++Y) { *Y = *X; }
                Y -= 2u*L;
                partial_sort_z(Y,L,k,ascend);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in partsort_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X1-=2u*L, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X1 -= 2u*L;
                    partial_sort_z(X1,L,k,ascend);
                    for (size_t l=L; l>0u; --l, ++X1, Y+=2u*K-1u) { *Y = *X1; *++Y = *++X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int partsort_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend, const size_t k)
{
    if (dim>3u) { fprintf(stderr,"error in partsort_inplace_s: dim must be in [0 3]\n"); return 1; }
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (k>L-1u) { fprintf(stderr,"error in partsort_inplace_s: k must be in [0 L-1]\n"); return 1; }
    
    // struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        partial_sort_s(X,L,k,ascend);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=L)
            {
                partial_sort_s(X,L,k,ascend);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in partsort_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X1, X+=K+1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    partial_sort_s(X1,L,k,ascend);
                    X1 += L-1u; X -= K;
                    for (size_t l=L; l>0u; --l, --X1, X-=K) { *X = *X1; }
                }
            }
            free(X1);
        }
    }

    // clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int partsort_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend, const size_t k)
{
    if (dim>3u) { fprintf(stderr,"error in partsort_inplace_d: dim must be in [0 3]\n"); return 1; }
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (k>L-1u) { fprintf(stderr,"error in partsort_inplace_d: k must be in [0 L-1]\n"); return 1; }

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        partial_sort_d(X,L,k,ascend);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=L)
            {
                partial_sort_d(X,L,k,ascend);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in partsort_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X1, X+=K+1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    partial_sort_d(X1,L,k,ascend);
                    X1 += L-1u; X -= K;
                    for (size_t l=L; l>0u; --l, --X1, X-=K) { *X = *X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int partsort_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend, const size_t k)
{
    if (dim>3u) { fprintf(stderr,"error in partsort_inplace_c: dim must be in [0 3]\n"); return 1; }
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (k>L-1u) { fprintf(stderr,"error in partsort_inplace_c: k must be in [0 L-1]\n"); return 1; }

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        partial_sort_c(X,L,k,ascend);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=2u*L)
            {
                partial_sort_c(X,L,k,ascend);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in partsort_inplace_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X1+=2, X+=2u*K+2u)
                {
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X1 -= 2u*L;
                    partial_sort_c(X1,L,k,ascend);
                    X1 += 2u*L-2u; X -= 2u*K;
                    for (size_t l=L; l>0u; --l, X1-=2, X-=2u*K) { *X = *X1; *(X+1) = *(X1+1); }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int partsort_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend, const size_t k)
{
    if (dim>3u) { fprintf(stderr,"error in partsort_inplace_z: dim must be in [0 3]\n"); return 1; }
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (k>L-1u) { fprintf(stderr,"error in partsort_inplace_z: k must be in [0 L-1]\n"); return 1; }

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        partial_sort_z(X,L,k,ascend);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=2u*L)
            {
                partial_sort_z(X,L,k,ascend);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in partsort_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X1+=2, X+=2u*K+2u)
                {
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X1 -= 2u*L;
                    partial_sort_z(X1,L,k,ascend);
                    X1 += 2u*L-2u; X -= 2u*K;
                    for (size_t l=L; l>0u; --l, X1-=2, X-=2u*K) { *X = *X1; *(X+1) = *(X1+1); }
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
