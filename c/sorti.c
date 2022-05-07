//Vec2vec operation.
//Sorts each vector in X along dim, and returns the indices.
//The indices are uints, but are returned as float or double.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "codee_math.h"
#include "cmpi_ascend.c"
#include "cmpi_descend.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int sorti_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sorti_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpi_ascend_s : cmpi_descend_s;

    FLT *XI;
    if (!(XI=(FLT *)malloc(L*sizeof(FLT)))) { fprintf(stderr,"error in sorti_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
        qsort(XI,L,sizeof(FLT),comp);
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
                qsort(XI,L,sizeof(FLT),comp);
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
                    qsort(XI,L,sizeof(FLT),comp);
                    for (size_t l=0u; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int sorti_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sorti_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpi_ascend_d : cmpi_descend_d;

    DBL *XI;
    if (!(XI=(DBL *)malloc(L*sizeof(DBL)))) { fprintf(stderr,"error in sorti_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
        qsort(XI,L,sizeof(DBL),comp);
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
                qsort(XI,L,sizeof(DBL),comp);
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
                    qsort(XI,L,sizeof(DBL),comp);
                    for (size_t l=0u; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int sorti_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sorti_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpi_ascend_c : cmpi_descend_c;

    CFLT *XI;
    if (!(XI=(CFLT *)malloc(L*sizeof(CFLT)))) { fprintf(stderr,"error in sorti_c: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = (float)l; }
        qsort(XI,L,sizeof(CFLT),comp);
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
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = (float)l; }
                qsort(XI,L,sizeof(CFLT),comp);
                for (size_t l=0u; l<L; ++l, ++Y) { *Y = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X, ++Y)
                {
                    for (size_t l=0u; l<L; ++l) { XI[l].r = X[2u*l*K]; XI[l].i = X[2u*l*K+1u]; XI[l].ind = (float)l; }
                    qsort(XI,L,sizeof(CFLT),comp);
                    for (size_t l=0u; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int sorti_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sorti_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpi_ascend_z : cmpi_descend_z;

    ZDBL *XI;
    if (!(XI=(ZDBL *)malloc(L*sizeof(ZDBL)))) { fprintf(stderr,"error in sorti_z: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = (double)l; }
        qsort(XI,L,sizeof(ZDBL),comp);
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
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = (double)l; }
                qsort(XI,L,sizeof(ZDBL),comp);
                for (size_t l=0u; l<L; ++l, ++Y) { *Y = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X, ++Y)
                {
                    for (size_t l=0u; l<L; ++l) { XI[l].r = X[2u*l*K]; XI[l].i = X[2u*l*K+1u]; XI[l].ind = (double)l; }
                    qsort(XI,L,sizeof(ZDBL),comp);
                    for (size_t l=0u; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int sorti_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sorti_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpi_ascend_s : cmpi_descend_s;

    FLT *XI;
    if (!(XI=(FLT *)malloc(L*sizeof(FLT)))) { fprintf(stderr,"error in sorti_inplace_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (float)l; }
        qsort(XI,L,sizeof(FLT),comp);
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
                qsort(XI,L,sizeof(FLT),comp);
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
                    qsort(XI,L,sizeof(FLT),comp);
                    for (size_t l=0u; l<L; ++l) { X[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int sorti_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in sorti_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpi_ascend_d : cmpi_descend_d;

    DBL *XI;
    if (!(XI=(DBL *)malloc(L*sizeof(DBL)))) { fprintf(stderr,"error in sorti_inplace_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (double)l; }
        qsort(XI,L,sizeof(DBL),comp);
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
                qsort(XI,L,sizeof(DBL),comp);
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
                    qsort(XI,L,sizeof(DBL),comp);
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
