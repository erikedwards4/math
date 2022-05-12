//Vec2vec operation.
//Sorts each vector in X along dim, and returns the indices.
//The indices are uints, but are returned as float or double.
//This uses the stdlib qsort function.

#include <stdio.h>
#include <stdlib.h>
#include "codee_math.h"
#include "cmpif.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int qsorti_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsorti_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpif_ascend_s : cmpif_descend_s;

    FLT_F *XI;
    if (!(XI=(FLT_F *)malloc(L*sizeof(FLT_F)))) { fprintf(stderr,"error in qsorti_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
        qsort(XI,L,sizeof(FLT_F),comp);
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
                qsort(XI,L,sizeof(FLT_F),comp);
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
                    qsort(XI,L,sizeof(FLT_F),comp);
                    for (size_t l=0u; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int qsorti_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsorti_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpif_ascend_d : cmpif_descend_d;

    DBL_D *XI;
    if (!(XI=(DBL_D *)malloc(L*sizeof(DBL_D)))) { fprintf(stderr,"error in qsorti_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
        qsort(XI,L,sizeof(DBL_D),comp);
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
                qsort(XI,L,sizeof(DBL_D),comp);
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
                    qsort(XI,L,sizeof(DBL_D),comp);
                    for (size_t l=0u; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int qsorti_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsorti_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpif_ascend_c : cmpif_descend_c;

    CFLT_F *XI;
    if (!(XI=(CFLT_F *)malloc(L*sizeof(CFLT_F)))) { fprintf(stderr,"error in qsorti_c: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = (float)l; }
        qsort(XI,L,sizeof(CFLT_F),comp);
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
                qsort(XI,L,sizeof(CFLT_F),comp);
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
                    qsort(XI,L,sizeof(CFLT_F),comp);
                    for (size_t l=0u; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int qsorti_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsorti_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpif_ascend_z : cmpif_descend_z;

    CDBL_D *XI;
    if (!(XI=(CDBL_D *)malloc(L*sizeof(CDBL_D)))) { fprintf(stderr,"error in qsorti_z: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = (double)l; }
        qsort(XI,L,sizeof(CDBL_D),comp);
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
                qsort(XI,L,sizeof(CDBL_D),comp);
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
                    qsort(XI,L,sizeof(CDBL_D),comp);
                    for (size_t l=0u; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int qsorti_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsorti_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpif_ascend_s : cmpif_descend_s;

    FLT_F *XI;
    if (!(XI=(FLT_F *)malloc(L*sizeof(FLT_F)))) { fprintf(stderr,"error in qsorti_inplace_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (float)l; }
        qsort(XI,L,sizeof(FLT_F),comp);
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
                qsort(XI,L,sizeof(FLT_F),comp);
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
                    qsort(XI,L,sizeof(FLT_F),comp);
                    for (size_t l=0u; l<L; ++l) { X[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int qsorti_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int ascend)
{
    if (dim>3u) { fprintf(stderr,"error in qsorti_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmpif_ascend_d : cmpif_descend_d;

    DBL_D *XI;
    if (!(XI=(DBL_D *)malloc(L*sizeof(DBL_D)))) { fprintf(stderr,"error in qsorti_inplace_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (double)l; }
        qsort(XI,L,sizeof(DBL_D),comp);
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
                qsort(XI,L,sizeof(DBL_D),comp);
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
                    qsort(XI,L,sizeof(DBL_D),comp);
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
