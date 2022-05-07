//Vec2scalar (reduction) operation.
//Does qselect algorithm to get kth largest element for each vector in X along dim.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

static size_t lomuto_partition_s (float *X, const size_t hi, size_t p, int largest);
static size_t lomuto_partition_d (double *X, const size_t hi, size_t p, int largest);
static size_t lomuto_partition_c (float *X, const size_t hi, size_t p, int largest);
static size_t lomuto_partition_z (double *X, const size_t hi, size_t p, int largest);
static float kselect_s (float *X, size_t hi, size_t k, int largest);
static double kselect_d (double *X, size_t hi, size_t k, int largest);
static size_t kselect_c (float *X, size_t hi, size_t k, int largest);
static size_t kselect_z (double *X, size_t hi, size_t k, int largest);


static size_t lomuto_partition_s (float *X, const size_t hi, size_t p, int largest)
{
    float x = X[p], pivot = x;
    X[p] = X[hi]; X[hi] = x;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i]>=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i]<=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    x = X[p]; X[p] = X[hi]; X[hi] = x;
    return p;
}


static size_t lomuto_partition_d (double *X, const size_t hi, size_t p, int largest)
{
    double x = X[p], pivot = x;
    X[p] = X[hi]; X[hi] = x;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i]>=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[i]<=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
        }
    }
    x = X[p]; X[p] = X[hi]; X[hi] = x;
    return p;
}


static size_t lomuto_partition_c (float *X, const size_t hi, size_t p, int largest)
{
    float xr = X[2u*p], xi = X[2u*p+1u], pivot = xr*xr + xi*xi;
    X[2u*p] = X[2u*hi]; X[2u*hi] = xr;
    X[2u*p+1u] = X[2u*hi+1u]; X[2u*hi+1u] = xi;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[2u*i]*X[2u*i]+X[2u*i+1u]*X[2u*i+1u]>=pivot+1e-5f)
            {
                xr = X[2u*i]; X[2u*i] = X[2u*p]; X[2u*p] = xr;
                xi = X[2u*i+1u]; X[2u*i+1u] = X[2u*p+1u]; X[2u*p+1u] = xi;
                ++p;
            }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[2u*i]*X[2u*i]+X[2u*i+1u]*X[2u*i+1u]<=pivot-1e-5f)
            {
                xr = X[2u*i]; X[2u*i] = X[2u*p]; X[2u*p] = xr;
                xi = X[2u*i+1u]; X[2u*i+1u] = X[2u*p+1u]; X[2u*p+1u] = xi;
                ++p;
            }
        }
    }
    xr = X[2u*p]; X[2u*p] = X[2u*hi]; X[2u*hi] = xr;
    xi = X[2u*p+1u]; X[2u*p+1u] = X[2u*hi+1u]; X[2u*hi+1u] = xi;
    return p;
}


static size_t lomuto_partition_z (double *X, const size_t hi, size_t p, int largest)
{
    double xr = X[2u*p], xi = X[2u*p+1u], pivot = xr*xr + xi*xi;
    X[2u*p] = X[2u*hi]; X[2u*hi] = xr;
    X[2u*p+1u] = X[2u*hi+1u]; X[2u*hi+1u] = xi;
    p = 0u;
    if (largest)
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[2u*i]*X[2u*i]+X[2u*i+1u]*X[2u*i+1u]>=pivot)
            {
                xr = X[2u*i]; X[2u*i] = X[2u*p]; X[2u*p] = xr;
                xi = X[2u*i+1u]; X[2u*i+1u] = X[2u*p+1u]; X[2u*p+1u] = xi;
                ++p;
            }
        }
    }
    else
    {
        for (size_t i=0u; i<hi; ++i)
        {
            if (X[2u*i]*X[2u*i]+X[2u*i+1u]*X[2u*i+1u]<=pivot)
            {
                xr = X[2u*i]; X[2u*i] = X[2u*p]; X[2u*p] = xr;
                xi = X[2u*i+1u]; X[2u*i+1u] = X[2u*p+1u]; X[2u*p+1u] = xi;
                ++p;
            }
        }
    }
    xr = X[2u*p]; X[2u*p] = X[2u*hi]; X[2u*hi] = xr;
    xi = X[2u*p+1u]; X[2u*p+1u] = X[2u*hi+1u]; X[2u*hi+1u] = xi;
    return p;
}


static float kselect_s (float *X, size_t hi, size_t k, int largest)
{
    while (1)
    {
        if (hi==0u) { return *X; }
        size_t p = hi/2u;
        p = lomuto_partition_s(X,hi,p,largest);
        if (k==p) { return *(X+k); }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
}


static double kselect_d (double *X, size_t hi, size_t k, int largest)
{
    while (1)
    {
        if (hi==0u) { return *X; }
        size_t p = hi/2u;
        p = lomuto_partition_d(X,hi,p,largest);
        if (k==p) { return *(X+k); }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
}


static size_t kselect_c (float *X, size_t hi, size_t k, int largest)
{
    size_t cnt = 0u;
    while (1)
    {
        if (hi==0u) { return cnt; }
        size_t p = hi/2u;
        p = lomuto_partition_c(X,hi,p,largest);
        if (k==p) { return cnt+2u*k; }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += 2u*p; cnt += 2u*p; hi -= p; k -= p; }
    }
}


static size_t kselect_z (double *X, size_t hi, size_t k, int largest)
{
    size_t cnt = 0u;
    while (1)
    {
        if (hi==0u) { return cnt; }
        size_t p = hi/2u;
        p = lomuto_partition_z(X,hi,p,largest);
        if (k==p) { return cnt+2u*k; }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += 2u*p; cnt += 2u*p; hi -= p; k -= p; }
    }
}


int qselect_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t k, const int largest)
{
    if (dim>3u) { fprintf(stderr,"error in qselect_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    // const char id = (largest) ? 'D' : 'I';
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (k>=L) { fprintf(stderr,"error in qselect_s: k must be in [0 L-1]\n"); return 1; }

    // struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        float *X1;
        if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in qselect_s: problem with malloc. "); perror("malloc"); return 1; }
        
        if (L==N)
        {
            for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
            X1 -= L;
            *Y = kselect_s(X1,L-1u,k,largest);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= L;
                    //if (LAPACKE_slasrt_work(id,(int)L,X1)) { fprintf(stderr,"error in qselect_s: problem with LAPACKE function\n"); }
                    //*Y = X1[k];
                    *Y = kselect_s(X1,L-1u,k,largest);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                    {
                        for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= L;
                        *Y = kselect_s(X1,L-1u,k,largest);
                    }
                }
            }
        }
        free(X1);
    }

    // clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int qselect_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t k, const int largest)
{
    if (dim>3u) { fprintf(stderr,"error in qselect_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (k>=L) { fprintf(stderr,"error in qselect_d: k must be in [0 L-1]\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        double *X1;
        if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in qselect_d: problem with malloc. "); perror("malloc"); return 1; }
        
        if (L==N)
        {
            for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
            X1 -= L;
            *Y = kselect_d(X1,L-1u,k,largest);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= L;
                    *Y = kselect_d(X1,L-1u,k,largest);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                    {
                        for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= L;
                        *Y = kselect_d(X1,L-1u,k,largest);
                    }
                }
            }
        }
        free(X1);
    }

    return 0;
}


int qselect_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t k, const int largest)
{
    if (dim>3u) { fprintf(stderr,"error in qselect_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (k>=L) { fprintf(stderr,"error in qselect_c: k must be in [0 L-1]\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        float *X1;
        if (!(X1=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in qselect_c: problem with malloc. "); perror("malloc"); return 1; }
        
        if (L==N)
        {
            for (size_t l=2u*L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
            X1 -= 2u*L;
            size_t i = kselect_c(X1,L-1u,k,largest);
            *Y = X1[i]; *(Y+1) = X1[i+1u];
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y+=2)
                {
                    for (size_t l=2u*L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2u*L;
                    size_t i = kselect_c(X1,L-1u,k,largest);
                    *Y = X1[i]; *(Y+1) = X1[i+1u];
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y+=2)
                    {
                        for (size_t l=L; l>0u; --l, X+=2u*K, X1+=2) { *X1 = *X; *(X1+1) = *(X+1); }
                        X1 -= 2u*L;
                        size_t i = kselect_c(X1,L-1u,k,largest);
                        *Y = X1[i]; *(Y+1) = X1[i+1u];
                    }
                }
            }
        }
        free(X1);
    }

    return 0;
}


int qselect_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t k, const int largest)
{
    if (dim>3u) { fprintf(stderr,"error in qselect_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (k>=L) { fprintf(stderr,"error in qselect_z: k must be in [0 L-1]\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        double *X1;
        if (!(X1=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in qselect_z: problem with malloc. "); perror("malloc"); return 1; }
        
        if (L==N)
        {
            for (size_t l=2u*L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
            X1 -= 2u*L;
            size_t i = kselect_z(X1,L-1u,k,largest);
            *Y = X1[i]; *(Y+1) = X1[i+1u];
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y+=2)
                {
                    for (size_t l=2u*L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2u*L;
                    size_t i = kselect_z(X1,L-1u,k,largest);
                    *Y = X1[i]; *(Y+1) = X1[i+1u];
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y+=2)
                    {
                        for (size_t l=L; l>0u; --l, X+=2u*K, X1+=2) { *X1 = *X; *(X1+1) = *(X+1); }
                        X1 -= 2u*L;
                        size_t i = kselect_z(X1,L-1u,k,largest);
                        *Y = X1[i]; *(Y+1) = X1[i+1u];
                    }
                }
            }
        }
        free(X1);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
