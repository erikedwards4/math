//Vec2scalar (reduction) operation.
//Gets index of median (50th percentile) for each vector in X along dim.
//This is the index with value closest to the 50th percentile,
//and rounds up for even-length vecs.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


typedef struct { float val; float ind; } FLT;
typedef struct { double val; double ind; } DBL;

int cmp_ascend_s (const void *a, const void *b);
int cmp_ascend_d (const void *a, const void *b);

int imed_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int imed_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


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


int imed_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in imed_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    FLT *XI;
    if (!(XI=(FLT *)malloc(L*sizeof(FLT)))) { fprintf(stderr,"error in imed_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0) {}
    else if (L==1)
    {
        float z = 0.0f;
        cblas_scopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
        qsort(XI,L,sizeof(FLT),cmp_ascend_s);
        *Y = XI[L/2].ind;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, ++Y)
            {
                for (size_t l=0; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
                qsort(XI,L,sizeof(FLT),cmp_ascend_s);
                *Y = XI[L/2].ind;
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    for (size_t l=0; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (float)l; }
                    qsort(XI,L,sizeof(FLT),cmp_ascend_s);
                    *Y = XI[L/2].ind;
                }
            }
        }
    }

    free(XI);
    return 0;
}


int imed_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in imed_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    DBL *XI;
    if (!(XI=(DBL *)malloc(L*sizeof(DBL)))) { fprintf(stderr,"error in imed_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        double z = 0.0;
        cblas_dcopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
        qsort(XI,L,sizeof(DBL),cmp_ascend_d);
        *Y = XI[L/2].ind;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, ++Y)
            {
                for (size_t l=0; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
                qsort(XI,L,sizeof(DBL),cmp_ascend_d);
                *Y = XI[L/2].ind;
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    for (size_t l=0; l<L; ++l) { XI[l].val = X[l]; XI[l].ind = (double)l; }
                    qsort(XI,L,sizeof(DBL),cmp_ascend_d);
                    *Y = XI[L/2].ind;
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
