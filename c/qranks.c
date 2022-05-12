//Vec2vec operation.
//Gets the rank (from 0 to L-1u) of the elements of each vector in X.
//Sorts each vector in X along dim, and returns the sorted rank of each element in X.
//This is the inverse of sorti, which returns the indices of sorted X to recover X,
//whereas this returns the indices of X to give sorted X.

#include <stdio.h>
#include <stdlib.h>
#include "codee_math.h"
#include "cmpi.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int qranks_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qranks_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpi_ascend_s : cmpi_descend_s;

    FLT_I *XI;
    if (!(XI=(FLT_I *)malloc(L*sizeof(FLT_I)))) { fprintf(stderr,"error in qranks_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = l; }
        qsort(XI,L,sizeof(FLT_I),comp);
        for (size_t l=0u; l<L; ++l) { Y[XI[l].ind] = (float)l; }
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
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = l; }
                qsort(XI,L,sizeof(FLT_I),comp);
                for (size_t l=0u; l<L; ++l) { Y[XI[l].ind] = (float)l; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X, ++Y)
                {
                    for (size_t l=0u; l<L; ++l) { XI[l].val = X[l*K]; XI[l].ind = l; }
                    qsort(XI,L,sizeof(FLT_I),comp);
                    for (size_t l=0u; l<L; ++l) { Y[K*XI[l].ind] = (float)l; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int qranks_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qranks_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpi_ascend_d : cmpi_descend_d;

    DBL_I *XI;
    if (!(XI=(DBL_I *)malloc(L*sizeof(DBL_I)))) { fprintf(stderr,"error in qranks_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = l; }
        qsort(XI,L,sizeof(DBL_I),comp);
        for (size_t l=0u; l<L; ++l) { Y[XI[l].ind] = (double)l; }
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
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = l; }
                qsort(XI,L,sizeof(DBL_I),comp);
                for (size_t l=0u; l<L; ++l) { Y[XI[l].ind] = (double)l; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X, ++Y)
                {
                    for (size_t l=0u; l<L; ++l) { XI[l].val = X[l*K]; XI[l].ind = l; }
                    qsort(XI,L,sizeof(DBL_I),comp);
                    for (size_t l=0u; l<L; ++l) { Y[K*XI[l].ind] = (double)l; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int qranks_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qranks_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpi_ascend_c : cmpi_descend_c;

    CFLT_I *XI;
    if (!(XI=(CFLT_I *)malloc(L*sizeof(CFLT_I)))) { fprintf(stderr,"error in qranks_c: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = l; }
        qsort(XI,L,sizeof(CFLT_I),comp);
        for (size_t l=0u; l<L; ++l) { Y[XI[l].ind] = (float)l; }
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
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = l; }
                qsort(XI,L,sizeof(CFLT_I),comp);
                for (size_t l=0u; l<L; ++l, ++Y) { Y[XI[l].ind] = (float)l; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X, ++Y)
                {
                    for (size_t l=0u; l<L; ++l) { XI[l].r = X[2u*l*K]; XI[l].i = X[2u*l*K+1u]; XI[l].ind = l; }
                    qsort(XI,L,sizeof(CFLT_I),comp);
                    for (size_t l=0u; l<L; ++l) { Y[K*XI[l].ind] = (float)l; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int qranks_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qranks_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpi_ascend_z : cmpi_descend_z;

    CDBL_I *XI;
    if (!(XI=(CDBL_I *)malloc(L*sizeof(CDBL_I)))) { fprintf(stderr,"error in qranks_z: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = l; }
        qsort(XI,L,sizeof(CDBL_I),comp);
        for (size_t l=0u; l<L; ++l) { Y[XI[l].ind] = (double)l; }
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
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = l; }
                qsort(XI,L,sizeof(CDBL_I),comp);
                for (size_t l=0u; l<L; ++l, ++Y) { Y[XI[l].ind] = (double)l; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X, ++Y)
                {
                    for (size_t l=0u; l<L; ++l) { XI[l].r = X[2u*l*K]; XI[l].i = X[2u*l*K+1u]; XI[l].ind = l; }
                    qsort(XI,L,sizeof(CDBL_I),comp);
                    for (size_t l=0u; l<L; ++l) { Y[K*XI[l].ind] = (double)l; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int qranks_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qranks_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpi_ascend_s : cmpi_descend_s;

    FLT_I *XI;
    if (!(XI=(FLT_I *)malloc(L*sizeof(FLT_I)))) { fprintf(stderr,"error in qranks_inplace_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = l; }
        qsort(XI,L,sizeof(FLT_I),comp);
        for (size_t l=0u; l<L; ++l) { X[XI[l].ind] = (float)l; }
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
                for (size_t l=0u; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = l; }
                qsort(XI,L,sizeof(FLT_I),comp);
                for (size_t l=0u; l<L; ++l) { X[XI[l].ind] = (float)l; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    for (size_t l=0u; l<L; ++l) { XI[l].val = X[l*K]; XI[l].ind = l; }
                    qsort(XI,L,sizeof(FLT_I),comp);
                    for (size_t l=0u; l<L; ++l) { X[K*XI[l].ind] = (float)l; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int qranks_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qranks_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpi_ascend_d : cmpi_descend_d;

    DBL_I *XI;
    if (!(XI=(DBL_I *)malloc(L*sizeof(DBL_I)))) { fprintf(stderr,"error in qranks_inplace_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = l; }
        qsort(XI,L,sizeof(DBL_I),comp);
        for (size_t l=0u; l<L; ++l) { X[XI[l].ind] = (double)l; }
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
                for (size_t l=0u; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = l; }
                qsort(XI,L,sizeof(DBL_I),comp);
                for (size_t l=0u; l<L; ++l) { X[XI[l].ind] = (double)l; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    for (size_t l=0u; l<L; ++l) { XI[l].val = X[l*K]; XI[l].ind = l; }
                    qsort(XI,L,sizeof(DBL_I),comp);
                    for (size_t l=0u; l<L; ++l) { X[K*XI[l].ind] = (double)l; }
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
