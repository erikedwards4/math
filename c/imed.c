//Vec2scalar (reduction) operation.
//Gets index of median (50th percentile) for each vector in X along dim.
//This is the index with value closest to the 50th percentile,
//and rounds up for even-length vecs.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "codee_math.h"
#include "cmpif.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int imed_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in imed_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    FLT_F *XI;
    if (!(XI=(FLT_F *)malloc(L*sizeof(FLT_F)))) { fprintf(stderr,"error in imed_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
        qsort(XI,L,sizeof(FLT_F),cmpif_ascend_s);
        *Y = XI[L/2u].ind;
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
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
                qsort(XI,L,sizeof(FLT_F),cmpif_ascend_s);
                *Y = XI[L/2u].ind;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K) { XI[l].val = *X; XI[l].ind = (float)l; }
                    qsort(XI,L,sizeof(FLT_F),cmpif_ascend_s);
                    *Y = XI[L/2u].ind;
                }
            }
        }
    }

    free(XI);
    return 0;
}


int imed_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in imed_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    DBL_D *XI;
    if (!(XI=(DBL_D *)malloc(L*sizeof(DBL_D)))) { fprintf(stderr,"error in imed_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
        qsort(XI,L,sizeof(DBL_D),cmpif_ascend_d);
        *Y = XI[L/2u].ind;
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
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
                qsort(XI,L,sizeof(DBL_D),cmpif_ascend_d);
                *Y = XI[L/2u].ind;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K) { XI[l].val = *X; XI[l].ind = (double)l; }
                    qsort(XI,L,sizeof(DBL_D),cmpif_ascend_d);
                    *Y = XI[L/2u].ind;
                }
            }
        }
    }

    free(XI);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
