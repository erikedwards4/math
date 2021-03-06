//Vec2scalar (reduction) operation.
//Gets index of pth percentile for each vector in X along dim.
//This is the index with value closest to the pth percentile,
//and rounds up for exact ties (0.5 case).


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


typedef struct { float val; float ind; } FLT;
typedef struct { double val; double ind; } DBL;

int cmp_ascend_s (const void *a, const void *b);
int cmp_ascend_d (const void *a, const void *b);

int iprctile_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p);
int iprctile_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p);


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


int iprctile_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const float p)
{
    if (dim>3) { fprintf(stderr,"error in iprctile_s: dim must be in [0 3]\n"); return 1; }
    if (p<0.0f || p>100.0f) { fprintf(stderr,"error in iprctile_s: p must be in [0 100]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    FLT *XI;
    if (!(XI=(FLT *)malloc(L*sizeof(FLT)))) { fprintf(stderr,"error in iprctile_s: problem with malloc. "); perror("malloc"); return 1; }

    //Get index closest to pth prctile after sorting
    const float p1 = (p/100.0f)*(L-1);
    const size_t i1 = (p<100.0f) ? (size_t)floorf(p1) : L-2;
    const float w2 = (p<100.0f) ? p1-floorf(p1) : 1.0f;
    const size_t ip = (w2<0.5f) ? i1 : i1+1;
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
        qsort(XI,L,sizeof(FLT),cmp_ascend_s);
        *Y = XI[ip].ind;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (float)l; }
                qsort(XI,L,sizeof(FLT),cmp_ascend_s);
                *Y = XI[ip].ind;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K) { XI[l].val = *X; XI[l].ind = (float)l; }
                    qsort(XI,L,sizeof(FLT),cmp_ascend_s);
                    *Y = XI[ip].ind;
                }
            }
        }
    }

    free(XI);
    return 0;
}


int iprctile_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const double p)
{
    if (dim>3) { fprintf(stderr,"error in iprctile_d: dim must be in [0 3]\n"); return 1; }
    if (p<0.0 || p>100.0) { fprintf(stderr,"error in iprctile_d: p must be in [0 100]"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    DBL *XI;
    if (!(XI=(DBL *)malloc(L*sizeof(DBL)))) { fprintf(stderr,"error in iprctile_d: problem with malloc. "); perror("malloc"); return 1; }

    //Get index closest to pth prctile after sorting
    const double p1 = (p/100.0)*(L-1);
    const size_t i1 = (p<100.0) ? (size_t)floor(p1) : L-2;
    const double w2 = (p<100.0) ? p1-floor(p1) : 1.0;
    const size_t ip = (w2<0.5) ? i1 : i1+1;
    
    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
        qsort(XI,L,sizeof(DBL),cmp_ascend_d);
        *Y = XI[ip].ind;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X) { XI[l].val = *X; XI[l].ind = (double)l; }
                qsort(XI,L,sizeof(DBL),cmp_ascend_d);
                *Y = XI[ip].ind;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X+=K) { XI[l].val = *X; XI[l].ind = (double)l; }
                    qsort(XI,L,sizeof(DBL),cmp_ascend_d);
                    *Y = XI[ip].ind;
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
