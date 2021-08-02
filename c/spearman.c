//Vecs2scalar operation for 2 inputs X1 and X2.
//Spearman rank correlation coefficient for each pair of vectors.

//This is quick implementation without fractional ranks and using d^2 formula (see Wikipedia).
//Thus, it gives correct answer if elements are unique within each of X1, X2.
//If elements are repeated within X1 or within X2, then answer is slightly off.
//I made the indices floats in case later want to use fractional ranks.

//Also maybe to do later: implement custom sort algorithm so don't need to copy into struct.

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

typedef struct { float val; size_t ind; } FLT;
typedef struct { double val; size_t ind; } DBL;

static int cmp_ascend_s (const void *a, const void *b);
static int cmp_ascend_d (const void *a, const void *b);

int spearman_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
int spearman_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);


static int cmp_ascend_s (const void *a, const void *b)
{
    const FLT x1 = *(const FLT *)a;
    const FLT x2 = *(const FLT *)b;
    if (x1.val>x2.val) { return 1; }
    else if (x2.val>x1.val) { return -1; }
    else { return 0; }
}


static int cmp_ascend_d (const void *a, const void *b)
{
	const DBL x1 = *(const DBL *)a;
    const DBL x2 = *(const DBL *)b;
    if (x1.val>x2.val) { return 1; }
    else if (x2.val>x1.val) { return -1; }
    else { return 0; }
}


int spearman_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in spearman_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in spearman_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    const int den = (int)(L*(L*L-1u));
    int d, dsm;
    int (*comp)(const void *, const void *) = cmp_ascend_s;

    int *r1, *r2;
    FLT *XI1, *XI2;
    if (!(r1=(int *)malloc(L*sizeof(int)))) { fprintf(stderr,"error in spearman_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(r2=(int *)malloc(L*sizeof(int)))) { fprintf(stderr,"error in spearman_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(XI1=(FLT *)malloc(L*sizeof(FLT)))) { fprintf(stderr,"error in spearman_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(XI2=(FLT *)malloc(L*sizeof(FLT)))) { fprintf(stderr,"error in spearman_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 1.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            XI1[l].val = *X1; XI1[l].ind = l;
            XI2[l].val = *X2; XI2[l].ind = l;
        }
        qsort(XI1,L,sizeof(FLT),comp); qsort(XI2,L,sizeof(FLT),comp);
        for (size_t l=0u; l<L; ++l) { r1[XI1[l].ind] = r2[XI2[l].ind] = (int)l; }
        dsm = 0;
        for (size_t l=0u; l<L; ++l, ++r1, ++r2) { d = *r1-*r2; dsm += d*d; }
        *Y = (float)(den-dsm*6) / (float)den;
        r1 -= L; r2 -= L;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? L : 0u, J2 = (L==N2) ? L : 0u;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2, r1-=L, r2-=L, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    XI1[l].val = *X1; XI1[l].ind = l;
                    XI2[l].val = *X2; XI2[l].ind = l;
                }
                qsort(XI1,L,sizeof(FLT),comp); qsort(XI2,L,sizeof(FLT),comp);
                for (size_t l=0u; l<L; ++l) { r1[XI1[l].ind] = r2[XI2[l].ind] = (int)l; }
                dsm = 0;
                for (size_t l=0u; l<L; ++l, ++r1, ++r2) { d = *r1-*r2; dsm += d*d; }
                *Y = (float)(den-dsm*6) / (float)den;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : B*(L-1u), I2 = (L==N2) ? 0u : B*(L-1u);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, r1-=L, r2-=L, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        XI1[l].val = *X1; XI1[l].ind = l;
                        XI2[l].val = *X2; XI2[l].ind = l;
                    }
                    qsort(XI1,L,sizeof(FLT),comp); qsort(XI2,L,sizeof(FLT),comp);
                    for (size_t l=0u; l<L; ++l) { r1[XI1[l].ind] = r2[XI2[l].ind] = (int)l; }
                    dsm = 0;
                    for (size_t l=0u; l<L; ++l, ++r1, ++r2) { d = *r1-*r2; dsm += d*d; }
                    *Y = (float)(den-dsm*6) / (float)den;
                }
            }
        }
    }

    free(XI1); free(XI2); free(r1); free(r2);
    return 0;
}


int spearman_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in spearman_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in spearman_d: vectors in X1 and X2 must have the same length\n"); return 1; }
    const int den = (int)(L*(L*L-1u));
    int d, dsm;
    int (*comp)(const void *, const void *) = cmp_ascend_d;

    int *r1, *r2;
    DBL *XI1, *XI2;
    if (!(r1=(int *)malloc(L*sizeof(int)))) { fprintf(stderr,"error in spearman_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(r2=(int *)malloc(L*sizeof(int)))) { fprintf(stderr,"error in spearman_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(XI1=(DBL *)malloc(L*sizeof(DBL)))) { fprintf(stderr,"error in spearman_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(XI2=(DBL *)malloc(L*sizeof(DBL)))) { fprintf(stderr,"error in spearman_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 1.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            XI1[l].val = *X1; XI1[l].ind = l;
            XI2[l].val = *X2; XI2[l].ind = l;
        }
        qsort(XI1,L,sizeof(DBL),comp); qsort(XI2,L,sizeof(DBL),comp);
        for (size_t l=0u; l<L; ++l) { r1[XI1[l].ind] = r2[XI2[l].ind] = (int)l; }
        dsm = 0;
        for (size_t l=0u; l<L; ++l, ++r1, ++r2) { d = *r1-*r2; dsm += d*d; }
        *Y = (double)(den-dsm*6) / (double)den;
        r1 -= L; r2 -= L;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? L : 0u, J2 = (L==N2) ? L : 0u;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2, r1-=L, r2-=L, ++Y)
            {
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    XI1[l].val = *X1; XI1[l].ind = l;
                    XI2[l].val = *X2; XI2[l].ind = l;
                }
                qsort(XI1,L,sizeof(DBL),comp); qsort(XI2,L,sizeof(DBL),comp);
                for (size_t l=0u; l<L; ++l) { r1[XI1[l].ind] = r2[XI2[l].ind] = (int)l; }
                dsm = 0;
                for (size_t l=0u; l<L; ++l, ++r1, ++r2) { d = *r1-*r2; dsm += d*d; }
                *Y = (double)(den-dsm*6) / (double)den;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : B*(L-1u), I2 = (L==N2) ? 0u : B*(L-1u);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, r1-=L, r2-=L, ++Y)
                {
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        XI1[l].val = *X1; XI1[l].ind = l;
                        XI2[l].val = *X2; XI2[l].ind = l;
                    }
                    qsort(XI1,L,sizeof(DBL),comp); qsort(XI2,L,sizeof(DBL),comp);
                    for (size_t l=0u; l<L; ++l) { r1[XI1[l].ind] = r2[XI2[l].ind] = (int)l; }
                    dsm = 0;
                    for (size_t l=0u; l<L; ++l, ++r1, ++r2) { d = *r1-*r2; dsm += d*d; }
                    *Y = (double)(den-dsm*6) / (double)den;
                }
            }
        }
    }

    free(XI1); free(XI2); free(r1); free(r2);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
