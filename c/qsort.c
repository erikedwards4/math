//Vec2vec operation.
//Sorts each vector in X along dim using quicksort.
//This is faster than insert_sort if the vectors are long (e.g., L>256).
//It is slightly faster than sort, which uses LAPACKE (same algorithm).
//This has in-place and not-in-place versions.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


#define HI_THRESH 128u
static size_t hoare_partition_s (float *X, const size_t N, const int ascend);
static size_t hoare_partition_d (double *X, const size_t N, const int ascend);
static size_t hoare_partition_c (float *X, const size_t N, const int ascend);
static size_t hoare_partition_z (double *X, const size_t N, const int ascend);
static void insertion_sort_s (float *X, const size_t N, const int ascend);
static void insertion_sort_d (double *X, const size_t N, const int ascend);
static void insertion_sort_c (float *X, const size_t N, const int ascend);
static void insertion_sort_z (double *X, const size_t N, const int ascend);
static void quicksort_s (float *X, const size_t N, const int ascend);
static void quicksort_d (double *X, const size_t N, const int ascend);
static void quicksort_c (float *X, const size_t N, const int ascend);
static void quicksort_z (double *X, const size_t N, const int ascend);


static size_t hoare_partition_s (float *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    float x, pivot = *X;
    if (ascend)
    {
        while (1)
        {
            while (X[i]<pivot) { ++i; }
            while (X[j]>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i]>pivot) { ++i; }
            while (X[j]<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


static size_t hoare_partition_d (double *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = N - 1u;
    double x, pivot = *X;
    if (ascend)
    {
        while (1)
        {
            while (X[i]<pivot) { ++i; }
            while (X[j]>pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
    else
    {
        while (1)
        {
            while (X[i]>pivot) { ++i; }
            while (X[j]<pivot) { --j; }
            if (j<=i) { return j+1u; }
            x = X[i]; X[i] = X[j]; X[j] = x;
            ++i; --j;
        }
    }
}


static size_t hoare_partition_c (float *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = 2u*N - 2u;
    float xr, xi, pivot = X[0]*X[0] + X[1]*X[1];
    if (ascend)
    {
        while (1)
        {
            while (X[i]*X[i]+X[i+1u]*X[i+1u]<pivot-1e-5f) { i += 2u; }
            while (X[j]*X[j]+X[j+1u]*X[j+1u]>pivot+1e-5f) { j -= 2u; }
            if (j<=i) { return j/2u+1u; }
            xr = X[i]; X[i] = X[j]; X[j] = xr;
            xi = X[i+1u]; X[i+1u] = X[j+1u]; X[j+1u] = xi;
            i += 2u; j -= 2u;
        }
    }
    else
    {
        while (1)
        {
            while (X[i]*X[i]+X[i+1u]*X[i+1u]>pivot+1e-5f) { i += 2u; }
            while (X[j]*X[j]+X[j+1u]*X[j+1u]<pivot-1e-5f) { j -= 2u; }
            if (j<=i) { return j/2u+1u; }
            xr = X[i]; X[i] = X[j]; X[j] = xr;
            xi = X[i+1u]; X[i+1u] = X[j+1u]; X[j+1u] = xi;
            i += 2u; j -= 2u;
        }
    }
}


static size_t hoare_partition_z (double *X, const size_t N, const int ascend)
{
    size_t i = 0u, j = 2u*N - 2u;
    double xr, xi, pivot = X[0]*X[0] + X[1]*X[1];
    if (ascend)
    {
        while (1)
        {
            while (X[i]*X[i]+X[i+1u]*X[i+1u]<pivot-1e-9) { i += 2u; }
            while (X[j]*X[j]+X[j+1u]*X[j+1u]>pivot+1e-9) { j -= 2u; }
            if (j<=i) { return j/2u+1u; }
            xr = X[i]; X[i] = X[j]; X[j] = xr;
            xi = X[i+1u]; X[i+1u] = X[j+1u]; X[j+1u] = xi;
            i += 2u; j -= 2u;
        }
    }
    else
    {
        while (1)
        {
            while (X[i]*X[i]+X[i+1u]*X[i+1u]>pivot+1e-9) { i += 2u; }
            while (X[j]*X[j]+X[j+1u]*X[j+1u]<pivot-1e-9) { j -= 2u; }
            if (j<=i) { return j/2u+1u; }
            xr = X[i]; X[i] = X[j]; X[j] = xr;
            xi = X[i+1u]; X[i+1u] = X[j+1u]; X[j+1u] = xi;
            i += 2u; j -= 2u;
        }
    }
}


static void insertion_sort_s(float *X, const size_t N, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<N; ++i)
        {
            size_t j = i; float x = *(X+1);
            while (j>0u && *X>x) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    else
    {
        for (size_t i=1u; i<N; ++i)
        {
            size_t j = i; float x = *(X+1);
            while (j>0u && *X<x) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    X -= N - 1u;
}


static void insertion_sort_d(double *X, const size_t N, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<N; ++i)
        {
            size_t j = i; double x = *(X+1);
            while (j>0u && *X>x) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    else
    {
        for (size_t i=1u; i<N; ++i)
        {
            size_t j = i; double x = *(X+1);
            while (j>0u && *X<x) { *(X+1) = *X; --X; --j; }
            *++X = x; X += i - j;
        }
    }
    X -= N - 1u;
}


static void insertion_sort_c(float *X, const size_t N, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<N; ++i)
        {
            size_t j = i; float xr = *(X+2), xi = *(X+3), x = xr*xr + xi*xi + 1e-5f;
            while (j>0u && *X**X+*(X+1)**(X+1)>x) { *(X+2) = *X; *(X+3) = *(X+1); X-=2; --j; }
            X += 2; *X = xr; *(X+1) = xi; X += 2u*(i-j);
        }
    }
    else
    {
        for (size_t i=1u; i<N; ++i)
        {
            size_t j = i; float xr = *(X+2), xi = *(X+3), x = xr*xr + xi*xi - 1e-5f;
            while (j>0u && *X**X+*(X+1)**(X+1)<x) { *(X+2) = *X; *(X+3) = *(X+1); X-=2; --j; }
            X += 2; *X = xr; *(X+1) = xi; X += 2u*(i-j);
        }
    }
    X -= 2u*N - 2u;
}


static void insertion_sort_z(double *X, const size_t N, const int ascend)
{
    if (ascend)
    {
        for (size_t i=1u; i<N; ++i)
        {
            size_t j = i; double xr = *(X+2), xi = *(X+3), x = xr*xr + xi*xi + 1e-9;
            while (j>0u && *X**X+*(X+1)**(X+1)>x) { *(X+2) = *X; *(X+3) = *(X+1); X-=2; --j; }
            X += 2; *X = xr; *(X+1) = xi; X += 2u*(i-j);
        }
    }
    else
    {
        for (size_t i=1u; i<N; ++i)
        {
            size_t j = i; double xr = *(X+2), xi = *(X+3), x = xr*xr + xi*xi - 1e-9;
            while (j>0u && *X**X+*(X+1)**(X+1)<x) { *(X+2) = *X; *(X+3) = *(X+1); X-=2; --j; }
            X += 2; *X = xr; *(X+1) = xi; X += 2u*(i-j);
        }
    }
    X -= 2u*N - 2u;
}


static void quicksort_s (float *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        // rearrange elements across pivot
        size_t p = hoare_partition_s(X, N, ascend);

        // recur on subarray containing elements that are < pivot
        if (p>HI_THRESH) { quicksort_s(X, p, ascend); }
        else { insertion_sort_s(X, p, ascend); }

        // recur on subarray containing elements that are > pivot
        if (N-p>HI_THRESH) { quicksort_s(X+p, N-p, ascend); }
        else { insertion_sort_s(X+p, N-p, ascend); }
    }
}


static void quicksort_d (double *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partition_d(X, N, ascend);

        if (p>HI_THRESH) { quicksort_d(X, p, ascend); }
        else { insertion_sort_d(X, p, ascend); }

        if (N-p>HI_THRESH) { quicksort_d(X+p, N-p, ascend); }
        else { insertion_sort_d(X+p, N-p, ascend); }
    }
}


static void quicksort_c (float *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partition_c(X, N, ascend);

        if (p>HI_THRESH) { quicksort_c(X, p, ascend); }
        else { insertion_sort_c(X, p, ascend); }

        if (N-p>HI_THRESH) { quicksort_c(X+2u*p, N-p, ascend); }
        else { insertion_sort_c(X+2u*p, N-p, ascend); }
    }
}


static void quicksort_z (double *X, const size_t N, const int ascend)
{
    if (N>1u)
    {
        size_t p = hoare_partition_z(X, N, ascend);

        if (p>HI_THRESH) { quicksort_z(X, p, ascend); }
        else { insertion_sort_z(X, p, ascend); }

        if (N-p>HI_THRESH) { quicksort_z(X+2u*p, N-p, ascend); }
        else { insertion_sort_z(X+2u*p, N-p, ascend); }
    }
}


int qsort_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsort_s: dim must be in [0 3]\n"); return 1; }

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
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in qsort_s: problem with malloc. "); perror("malloc"); return 1; }
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


int qsort_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsort_d: dim must be in [0 3]\n"); return 1; }

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
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in qsort_d: problem with malloc. "); perror("malloc"); return 1; }
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


int qsort_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsort_c: dim must be in [0 3]\n"); return 1; }

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
            if (!(X1=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in qsort_c: problem with malloc. "); perror("malloc"); return 1; }
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


int qsort_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsort_z: dim must be in [0 3]\n"); return 1; }

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
            if (!(X1=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in qsort_z: problem with malloc. "); perror("malloc"); return 1; }
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


int qsort_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsort_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    
    struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

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
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in qsort_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
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

    clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int qsort_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsort_inplace_d: dim must be in [0 3]\n"); return 1; }

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
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in qsort_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
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


int qsort_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsort_inplace_c: dim must be in [0 3]\n"); return 1; }

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
            if (!(X1=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in qsort_inplace_c: problem with malloc. "); perror("malloc"); return 1; }
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


int qsort_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsort_inplace_z: dim must be in [0 3]\n"); return 1; }

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
            if (!(X1=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in qsort_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
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
