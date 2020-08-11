//Vec2vec operation.
//Sorts each vector in X along dim, and returns the indices.
//The indices are uints, but are returned as float or double.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

//It is faster to keep the inds as the same data type as the vals.
typedef struct { float val; float ind; } FLT;
typedef struct { double val; double ind; } DBL;
typedef struct { float r; float i; float ind; } CFLT;
typedef struct { double r; double i; double ind; } ZDBL;

int cmp_ascend_s (const void *a, const void *b);
int cmp_ascend_d (const void *a, const void *b);
int cmp_ascend_c (const void *a, const void *b);
int cmp_ascend_z (const void *a, const void *b);

int cmp_descend_s (const void *a, const void *b);
int cmp_descend_d (const void *a, const void *b);
int cmp_descend_c (const void *a, const void *b);
int cmp_descend_z (const void *a, const void *b);

int sorti_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sorti_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sorti_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sorti_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);

int sorti_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sorti_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);


int cmp_ascend_s (const void *a, const void *b)
{
    const FLT x1 = *(const FLT *)a;
    const FLT x2 = *(const FLT *)b;
    if (x1.val>x2.val) { return 1; }
    else if (x2.val>x1.val) { return -1; }
    else { return 0; }
}

int cmp_ascend_d (const void *a, const void *b)
{
	const DBL x1 = *(const DBL *)a;
    const DBL x2 = *(const DBL *)b;
    if (x1.val>x2.val) { return 1; }
    else if (x2.val>x1.val) { return -1; }
    else { return 0; }
}

int cmp_ascend_c (const void *a, const void *b)
{
    const CFLT x1 = *(const CFLT *)a;
    const CFLT x2 = *(const CFLT *)b;
	const float sq1 = x1.r*x1.r + x1.i*x1.i;
    const float sq2 = x2.r*x2.r + x2.i*x2.i;
	if (sq1>sq2) { return 1; }
    else if (sq2>sq1) { return -1; }
	else if (atan2f(x1.i,x1.r) > atan2f(x2.i,x2.r)) { return 1; }
    else if (atan2f(x2.i,x2.r) > atan2f(x1.i,x1.r)) { return -1; }
    else { return 0; }
}

int cmp_ascend_z (const void *a, const void *b)
{
	const ZDBL x1 = *(const ZDBL *)a;
    const ZDBL x2 = *(const ZDBL *)b;
	const double sq1 = x1.r*x1.r + x1.i*x1.i;
    const double sq2 = x2.r*x2.r + x2.i*x2.i;
	if (sq1>sq2) { return 1; }
    else if (sq2>sq1) { return -1; }
	else if (atan2(x1.i,x1.r) > atan2(x2.i,x2.r)) { return 1; }
    else if (atan2(x2.i,x2.r) > atan2(x1.i,x1.r)) { return -1; }
    else { return 0; }
}


int cmp_descend_s (const void *a, const void *b)
{
    const FLT x1 = *(const FLT *)a;
    const FLT x2 = *(const FLT *)b;
    if (x1.val>x2.val) { return -1; }
    else if (x2.val>x1.val) { return 1; }
    else { return 0; }
}

int cmp_descend_d (const void *a, const void *b)
{
	const DBL x1 = *(const DBL *)a;
    const DBL x2 = *(const DBL *)b;
    if (x1.val>x2.val) { return -1; }
    else if (x2.val>x1.val) { return 1; }
    else { return 0; }
}

int cmp_descend_c (const void *a, const void *b)
{
    const CFLT x1 = *(const CFLT *)a;
    const CFLT x2 = *(const CFLT *)b;
	const float sq1 = x1.r*x1.r + x1.i*x1.i;
    const float sq2 = x2.r*x2.r + x2.i*x2.i;
	if (sq1>sq2) { return -1; }
    else if (sq2>sq1) { return 1; }
	else if (atan2f(x1.i,x1.r) > atan2f(x2.i,x2.r)) { return -1; }
    else if (atan2f(x2.i,x2.r) > atan2f(x1.i,x1.r)) { return 1; }
    else { return 0; }
}

int cmp_descend_z (const void *a, const void *b)
{
	const ZDBL x1 = *(const ZDBL *)a;
    const ZDBL x2 = *(const ZDBL *)b;
	const double sq1 = x1.r*x1.r + x1.i*x1.i;
    const double sq2 = x2.r*x2.r + x2.i*x2.i;
	if (sq1>sq2) { return -1; }
    else if (sq2>sq1) { return 1; }
	else if (atan2(x1.i,x1.r) > atan2(x2.i,x2.r)) { return -1; }
    else if (atan2(x2.i,x2.r) > atan2(x1.i,x1.r)) { return 1; }
    else { return 0; }
}


int sorti_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sorti_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_s : cmp_descend_s;

    FLT *XI;
    if (!(XI=(FLT *)malloc(L*sizeof(FLT)))) { fprintf(stderr,"error in sorti_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
        qsort(XI,L,sizeof(FLT),comp);
        for (size_t l=0; l<L; ++l, ++Y) { *Y = XI[l].ind; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t l=0; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
                qsort(XI,L,sizeof(FLT),comp);
                for (size_t l=0; l<L; ++l, ++Y) { *Y = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    for (size_t l=0; l<L; ++l) { XI[l].val = X[l*K]; XI[l].ind = (float)l; }
                    qsort(XI,L,sizeof(FLT),comp);
                    for (size_t l=0; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int sorti_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sorti_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_d : cmp_descend_d;

    DBL *XI;
    if (!(XI=(DBL *)malloc(L*sizeof(DBL)))) { fprintf(stderr,"error in sorti_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
        qsort(XI,L,sizeof(DBL),comp);
        for (size_t l=0; l<L; ++l, ++Y) { *Y = XI[l].ind; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t l=0; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
                qsort(XI,L,sizeof(DBL),comp);
                for (size_t l=0; l<L; ++l, ++Y) { *Y = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    for (size_t l=0; l<L; ++l) { XI[l].val = X[l*K]; XI[l].ind = (double)l; }
                    qsort(XI,L,sizeof(DBL),comp);
                    for (size_t l=0; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int sorti_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sorti_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_c : cmp_descend_c;

    CFLT *XI;
    if (!(XI=(CFLT *)malloc(L*sizeof(CFLT)))) { fprintf(stderr,"error in sorti_c: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = (float)l; }
        qsort(XI,L,sizeof(CFLT),comp);
        for (size_t l=0; l<L; ++l, ++Y) { *Y = XI[l].ind; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t l=0; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = (float)l; }
                qsort(XI,L,sizeof(CFLT),comp);
                for (size_t l=0; l<L; ++l, ++Y) { *Y = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    for (size_t l=0; l<L; ++l) { XI[l].r = X[2*l*K]; XI[l].i = X[2*l*K+1]; XI[l].ind = (float)l; }
                    qsort(XI,L,sizeof(CFLT),comp);
                    for (size_t l=0; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int sorti_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sorti_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_z : cmp_descend_z;

    ZDBL *XI;
    if (!(XI=(ZDBL *)malloc(L*sizeof(ZDBL)))) { fprintf(stderr,"error in sorti_z: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = (double)l; }
        qsort(XI,L,sizeof(ZDBL),comp);
        for (size_t l=0; l<L; ++l, ++Y) { *Y = XI[l].ind; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t l=0; l<L; ++l, ++X) { XI[l].r = *X; XI[l].i = *++X; XI[l].ind = (double)l; }
                qsort(XI,L,sizeof(ZDBL),comp);
                for (size_t l=0; l<L; ++l, ++Y) { *Y = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    for (size_t l=0; l<L; ++l) { XI[l].r = X[2*l*K]; XI[l].i = X[2*l*K+1]; XI[l].ind = (double)l; }
                    qsort(XI,L,sizeof(ZDBL),comp);
                    for (size_t l=0; l<L; ++l) { Y[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int sorti_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sorti_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_s : cmp_descend_s;

    FLT *XI;
    if (!(XI=(FLT *)malloc(L*sizeof(FLT)))) { fprintf(stderr,"error in sorti_inplace_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (float)l; }
        qsort(XI,L,sizeof(FLT),comp);
        for (size_t l=0; l<L; ++l, ++X) { *X = XI[l].ind; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t l=0; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (float)l; }
                qsort(XI,L,sizeof(FLT),comp);
                for (size_t l=0; l<L; ++l, ++X) { *X = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X)
                {
                    for (size_t l=0; l<L; ++l) { XI[l].val = X[l*K]; XI[l].ind = (float)l; }
                    qsort(XI,L,sizeof(FLT),comp);
                    for (size_t l=0; l<L; ++l) { X[l*K] = XI[l].ind; }
                }
            }
        }
    }

    free(XI);
    return 0;
}


int sorti_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sorti_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_d : cmp_descend_d;

    DBL *XI;
    if (!(XI=(DBL *)malloc(L*sizeof(DBL)))) { fprintf(stderr,"error in sorti_inplace_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (double)l; }
        qsort(XI,L,sizeof(DBL),comp);
        for (size_t l=0; l<L; ++l, ++X) { *X = XI[l].ind; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t l=0; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (double)l; }
                qsort(XI,L,sizeof(DBL),comp);
                for (size_t l=0; l<L; ++l, ++X) { *X = XI[l].ind; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X)
                {
                    for (size_t l=0; l<L; ++l) { XI[l].val = X[l*K]; XI[l].ind = (double)l; }
                    qsort(XI,L,sizeof(DBL),comp);
                    for (size_t l=0; l<L; ++l) { X[l*K] = XI[l].ind; }
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
