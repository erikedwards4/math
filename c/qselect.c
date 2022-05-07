//Vec2scalar (reduction) operation.
//Does qselect algorithm to get kth largest element for each vector in X along dim.

#include <stdio.h>
#include <lapacke.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

static size_t lomuto_partition_s (float *X, const size_t hi, size_t p);
static float kselect_s (float *X, size_t hi, size_t k);
// static float kselect_s (float *X, const size_t lo, const size_t hi, const size_t k);


static size_t lomuto_partition_s (float *X, const size_t hi, size_t p)
{
    float x = X[p], pivot = x;
    X[p] = X[hi]; X[hi] = x;
    p = 0u;
    for (size_t i=0u; i<hi; ++i)
    {
        if (X[i]<=pivot) { x = X[i]; X[i] = X[p]; X[p] = x; ++p; }
    }
    x = X[p]; X[p] = X[hi]; X[hi] = x;
    return p;
}


static float kselect_s (float *X, size_t hi, size_t k)
{
    while (1)
    {
        if (hi==0u) { return *X; }
        size_t p = hi/2u;
        p = lomuto_partition_s(X,hi,p);
        if (k==p) { return *(X+k); }
        else if (k<p) { hi = p - 1u; }
        else { ++p; X += p; hi -= p; k -= p; }
    }
}

// static float kselect_s (float *X, const size_t lo, const size_t hi, const size_t k)
// {
//     if (hi==lo) { return X[lo]; }
//     size_t p = lo + (hi-lo)/2u;  //Init p (pivot index) between lo and hi
//     p = lomuto_partition_s(X,hi,p);
//     if (k==p) { return X[k]; }
//     else if (k<p) { return kselect_s(X,lo,p-1u,k); }
//     else { return kselect_s(X,p+1u,hi,k); }
// }


int qselect_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t k, const int largest)
{
    if (dim>3u) { fprintf(stderr,"error in qselect_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    // const char id = (largest) ? 'D' : 'I';
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (k>=L) { fprintf(stderr,"error in qselect_s: k must be in [1 L]\n"); return 1; }

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
            //if (LAPACKE_slasrt_work(id,(int)L,X1)) { fprintf(stderr,"error in qselect_s: problem with LAPACKE function\n"); }
            //*Y = X1[k];
            *Y = kselect_s(X1,L-1u,k);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= L;
                    //if (LAPACKE_slasrt_work(id,(int)L,X1)) { fprintf(stderr,"error in qselect_s: problem with LAPACKE function\n"); }
                    //*Y = X1[k];
                    *Y = kselect_s(X1,L-1u,k);
                }
                clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                    {
                        for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= L;
                        //if (LAPACKE_slasrt_work(id,(int)L,X1)) { fprintf(stderr,"error in qselect_s: problem with LAPACKE function\n"); }
                        //*Y = X1[k];
                        *Y = kselect_s(X1,L-1u,k);
                    }
                }
            }
        }
        free(X1);
    }

    return 0;
}


int qselect_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t k, const int largest)
{
    if (dim>3u) { fprintf(stderr,"error in qselect_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const char id = (largest) ? 'D' : 'I';

    if (k<1u || k>L) { fprintf(stderr,"error in qselect_d: k must be in [1 L]\n"); return 1; }

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
            if (LAPACKE_dlasrt_work(id,(int)L,X1)) { fprintf(stderr,"error in qselect_d: problem with LAPACKE function\n"); }
            *Y = X1[k-1u];
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
                    if (LAPACKE_dlasrt_work(id,(int)L,X1)) { fprintf(stderr,"error in qselect_d: problem with LAPACKE function\n"); }
                    *Y = X1[k-1u];
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
                        if (LAPACKE_dlasrt_work(id,(int)L,X1)) { fprintf(stderr,"error in qselect_d: problem with LAPACKE function\n"); }
                        *Y = X1[k-1u];
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
