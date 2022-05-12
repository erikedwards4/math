//Vec2vec operation.
//Sorts each vector in X along dim, and returns the indices.
//The indices are uints, but are returned as float or double.

#include <stdio.h>
#include <stdlib.h>
#include "codee_math.h"
#include "insertion_sortif.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int insert_sorti_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in insert_sorti_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    FLT_F *XI;
    if (!(XI=(FLT_F *)malloc(L*sizeof(FLT_F)))) { fprintf(stderr,"error in insert_sorti_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
        insertion_sortif_s(XI,L,ascend);
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = XI[l].ind; }
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
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
                insertion_sortif_s(XI,L,ascend);
                for (size_t l=0u; l<L; ++l, ++Y) { *Y = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X, ++Y)
                {
                    for (size_t l=0u; l<L; ++l) { XI[l].val = X[l*K]; XI[l].ind = (float)l; }
                    insertion_sortif_s(XI,L,ascend);
                    for (size_t l=0u; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int insert_sorti_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in insert_sorti_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    DBL_D *XI;
    if (!(XI=(DBL_D *)malloc(L*sizeof(DBL_D)))) { fprintf(stderr,"error in insert_sorti_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
        insertion_sortif_d(XI,L,ascend);
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = XI[l].ind; }
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
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
                insertion_sortif_d(XI,L,ascend);
                for (size_t l=0u; l<L; ++l, ++Y) { *Y = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X, ++Y)
                {
                    for (size_t l=0u; l<L; ++l) { XI[l].val = X[l*K]; XI[l].ind = (double)l; }
                    insertion_sortif_d(XI,L,ascend);
                    for (size_t l=0u; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int insert_sorti_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in insert_sorti_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    FLT_F *XI;
    if (!(XI=(FLT_F *)malloc(L*sizeof(FLT_F)))) { fprintf(stderr,"error in insert_sorti_c: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, X+=2) { XI[l].val = *X**X+*(X+1)**(X+1); XI[l].ind = (float)l; }
        insertion_sortif_s(XI,L,ascend);
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = XI[l].ind; }
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
                for (size_t l=0u; l<L; ++l, X+=2) { XI[l].val = *X**X+*(X+1)**(X+1); XI[l].ind = (float)l; }
                insertion_sortif_s(XI,L,ascend);
                for (size_t l=0u; l<L; ++l, ++Y) { *Y = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<L; ++l, X+=2u*K) { XI[l].val = *X**X+*(X+1)**(X+1); XI[l].ind = (float)l; }
                    insertion_sortif_s(XI,L,ascend);
                    for (size_t l=0u; l<L; ++l, Y+=K) { *Y = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int insert_sorti_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in insert_sorti_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    DBL_D *XI;
    if (!(XI=(DBL_D *)malloc(L*sizeof(DBL_D)))) { fprintf(stderr,"error in insert_sorti_z: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, X+=2) { XI[l].val = *X**X+*(X+1)**(X+1); XI[l].ind = (double)l; }
        insertion_sortif_d(XI,L,ascend);
        for (size_t l=0u; l<L; ++l, ++Y) { *Y = XI[l].ind; }
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
                for (size_t l=0u; l<L; ++l, X+=2) { XI[l].val = *X**X+*(X+1)**(X+1); XI[l].ind = (double)l; }
                insertion_sortif_d(XI,L,ascend);
                for (size_t l=0u; l<L; ++l, ++Y) { *Y = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<L; ++l, X+=2u*K) { XI[l].val = *X**X+*(X+1)**(X+1); XI[l].ind = (double)l; }
                    insertion_sortif_d(XI,L,ascend);
                    for (size_t l=0u; l<L; ++l, Y+=K) { *Y = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int insert_sorti_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in insert_sorti_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    FLT_F *XI;
    if (!(XI=(FLT_F *)malloc(L*sizeof(FLT_F)))) { fprintf(stderr,"error in insert_sorti_inplace_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (float)l; }
        insertion_sortif_s(XI,L,ascend);
        for (size_t l=0u; l<L; ++l, ++X) { *X = XI[l].ind; }
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
                for (size_t l=0u; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (float)l; }
                insertion_sortif_s(XI,L,ascend);
                for (size_t l=0u; l<L; ++l, ++X) { *X = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    for (size_t l=0u; l<L; ++l) { XI[l].val = X[l*K]; XI[l].ind = (float)l; }
                    insertion_sortif_s(XI,L,ascend);
                    for (size_t l=0u; l<L; ++l) { X[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int insert_sorti_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in insert_sorti_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    DBL_D *XI;
    if (!(XI=(DBL_D *)malloc(L*sizeof(DBL_D)))) { fprintf(stderr,"error in insert_sorti_inplace_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (double)l; }
        insertion_sortif_d(XI,L,ascend);
        for (size_t l=0u; l<L; ++l, ++X) { *X = XI[l].ind; }
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
                for (size_t l=0u; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (double)l; }
                insertion_sortif_d(XI,L,ascend);
                for (size_t l=0u; l<L; ++l, ++X) { *X = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    for (size_t l=0u; l<L; ++l) { XI[l].val = X[l*K]; XI[l].ind = (double)l; }
                    insertion_sortif_d(XI,L,ascend);
                    for (size_t l=0u; l<L; ++l) { X[l*K] = XI[l].ind; }
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
