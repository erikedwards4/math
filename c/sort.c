//Vec2vec operation.
//Sorts each vector in X along dim using quicksort.
//This is faster than insert_sort if the vectors are long (e.g., L>256).
//It is slightly faster than lapacke_sort, which uses LAPACKE (same algorithm).
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <stdlib.h>
#include "codee_math.h"
#include "quicksort.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int sort_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        Y -= L;
        quicksort_s(Y,L,ascend);
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
                quicksort_s(Y,L,ascend);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in sort_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, X1-=L, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    quicksort_s(X1,L,ascend);
                    for (size_t l=L; l>0u; --l, ++X1, Y+=K) { *Y = *X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        quicksort_d(Y,L,ascend);
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
                quicksort_d(Y,L,ascend);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in sort_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, X1-=L, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    quicksort_d(X1,L,ascend);
                    for (size_t l=L; l>0u; --l, ++X1, Y+=K) { *Y = *X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=2u*L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        quicksort_c(Y,L,ascend);
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
                quicksort_c(Y,L,ascend);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in sort_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X1-=2u*L, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X1 -= 2u*L;
                    quicksort_c(X1,L,ascend);
                    for (size_t l=L; l>0u; --l, ++X1, Y+=2u*K-1u) { *Y = *X1; *++Y = *++X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=2u*L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        quicksort_z(Y,L,ascend);
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
                quicksort_z(Y,L,ascend);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in sort_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X1-=2u*L, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X1 -= 2u*L;
                    quicksort_z(X1,L,ascend);
                    for (size_t l=L; l>0u; --l, ++X1, Y+=2u*K-1u) { *Y = *X1; *++Y = *++X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    
    // struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        quicksort_s(X,L,ascend);
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
                quicksort_s(X,L,ascend);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in sort_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X1, X+=K+1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    quicksort_s(X1,L,ascend);
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


int sort_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        quicksort_d(X,L,ascend);
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
                quicksort_d(X,L,ascend);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in sort_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X1, X+=K+1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    quicksort_d(X1,L,ascend);
                    X1 += L-1u; X -= K;
                    for (size_t l=L; l>0u; --l, --X1, X-=K) { *X = *X1; }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_inplace_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        quicksort_c(X,L,ascend);
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
                quicksort_c(X,L,ascend);
            }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in sort_inplace_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X1+=2, X+=2u*K+2u)
                {
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X1 -= 2u*L;
                    quicksort_c(X1,L,ascend);
                    X1 += 2u*L-2u; X -= 2u*K;
                    for (size_t l=L; l>0u; --l, X1-=2, X-=2u*K) { *X = *X1; *(X+1) = *(X1+1); }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int sort_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sort_inplace_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        quicksort_z(X,L,ascend);
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
                quicksort_z(X,L,ascend);
            }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in sort_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X1+=2, X+=2u*K+2u)
                {
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X1 -= 2u*L;
                    quicksort_z(X1,L,ascend);
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
