//Vecs2scalar operation for 2 inputs X1 and X2.
//Spearman rank correlation coefficient for each pair of vectors.

//This is quick implementation without fractional ranks and using d^2 formula (see Wikipedia).
//Thus, it gives correct answer if elements are unique within each of X1, X2.
//If elements are repeated within X1 or within X2, then answer is slightly off.
//I made the indices floats in case later want to use fractional ranks.

//Also maybe to do later: implement custom sort algorithm so don't need to copy into struct.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <lapacke.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

typedef struct { float val; size_t ind; } FLT;
typedef struct { double val; size_t ind; } DBL;

int cmp_ascend_s (const void *a, const void *b);
int cmp_ascend_d (const void *a, const void *b);

int spearman_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int spearman_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);


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


int spearman_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in spearman_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in spearman_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float den = 6.0f/(L*(L*L-1));
    float d, dsm;
    int (*comp)(const void *, const void *) = cmp_ascend_s;

    float *r1, *r2;
    FLT *XI1, *XI2;
    if (!(r1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in spearman_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(r2=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in spearman_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(XI1=(FLT *)malloc(L*sizeof(FLT)))) { fprintf(stderr,"error in spearman_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(XI2=(FLT *)malloc(L*sizeof(FLT)))) { fprintf(stderr,"error in spearman_s: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        const float o = 1.0f;
        cblas_scopy((int)N,&o,0,Y,1);
    }
    else if (L==N)
    {
        //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
        //cblas_scopy((int)L,X1,1,(float *)XI1,2); //this works, but not faster
        //cblas_scopy((int)L,X2,1,(float *)XI2,2);
        for (size_t l=0; l<L; ++l)
        {
            XI1[l].val = X1[l]; XI1[l].ind = l;
            XI2[l].val = X2[l]; XI2[l].ind = l;
        }
        qsort(XI1,L,sizeof(FLT),comp); qsort(XI2,L,sizeof(FLT),comp);
        for (size_t l=0; l<L; ++l) { r1[XI1[l].ind] = r2[XI2[l].ind] = (float)l; }
        dsm = 0.0f; for (size_t l=0; l<L; ++l) { d = r1[l]-r2[l]; dsm += d*d; }
        *Y = 1.0f - dsm*den;
        //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? L : 0, J2 = (L==N2) ? L : 0;
            for (size_t v=0; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                for (size_t l=0; l<L; ++l)
                {
                    XI1[l].val = *X1++; XI1[l].ind = l;
                    XI2[l].val = *X2++; XI2[l].ind = l;
                }
                qsort(XI1,L,sizeof(FLT),comp); qsort(XI2,L,sizeof(FLT),comp);
                for (size_t l=0; l<L; ++l) { r1[XI1[l].ind] = r2[XI2[l].ind] = (float)l; }
                dsm = 0.0f; for (size_t l=0; l<L; ++l) { d = r1[l]-r2[l]; dsm += d*d; }
                *Y = 1.0f - dsm*den;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            const size_t K1 = (L==N1) ? 1 : K, K2 = (L==N2) ? 1 : K;
            const size_t I1 = (L==N1) ? 0 : B*(L-1), I2 = (L==N2) ? 0 : B*(L-1);
            for (size_t g=0; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    for (size_t l=0; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        XI1[l].val = *X1; XI1[l].ind = l;
                        XI2[l].val = *X2; XI2[l].ind = l;
                    }
                    qsort(XI1,L,sizeof(FLT),comp); qsort(XI2,L,sizeof(FLT),comp);
                    for (size_t l=0; l<L; ++l) { r1[XI1[l].ind] = r2[XI2[l].ind] = (float)l; }
                    dsm = 0.0f; for (size_t l=0; l<L; ++l) { d = r1[l]-r2[l]; dsm += d*d; }
                    *Y = 1.0f - dsm*den;
                }
            }
        }
    }

    free(XI1); free(XI2); free(r1); free(r2);
    return 0;
}


int spearman_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in spearman_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in spearman_d: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double den = 6.0/(L*(L*L-1));
    double d, dsm;
    int (*comp)(const void *, const void *) = cmp_ascend_d;

    double *r1, *r2;
    DBL *XI1, *XI2;
    if (!(r1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in spearman_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(r2=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in spearman_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(XI1=(DBL *)malloc(L*sizeof(DBL)))) { fprintf(stderr,"error in spearman_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(XI2=(DBL *)malloc(L*sizeof(DBL)))) { fprintf(stderr,"error in spearman_d: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        const double o = 1.0;
        cblas_dcopy((int)N,&o,0,Y,1);
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l)
        {
            XI1[l].val = X1[l]; XI1[l].ind = l;
            XI2[l].val = X2[l]; XI2[l].ind = l;
        }
        qsort(XI1,L,sizeof(DBL),comp); qsort(XI2,L,sizeof(DBL),comp);
        for (size_t l=0; l<L; ++l) { r1[XI1[l].ind] = r2[XI2[l].ind] = (double)l; }
        dsm = 0.0; for (size_t l=0; l<L; ++l) { d = r1[l]-r2[l]; dsm += d*d; }
        *Y = 1.0 - dsm*den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? L : 0, J2 = (L==N2) ? L : 0;
            for (size_t v=0; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                for (size_t l=0; l<L; ++l)
                {
                    XI1[l].val = *X1++; XI1[l].ind = l;
                    XI2[l].val = *X2++; XI2[l].ind = l;
                }
                qsort(XI1,L,sizeof(DBL),comp); qsort(XI2,L,sizeof(DBL),comp);
                for (size_t l=0; l<L; ++l) { r1[XI1[l].ind] = r2[XI2[l].ind] = (double)l; }
                dsm = 0.0; for (size_t l=0; l<L; ++l) { d = r1[l]-r2[l]; dsm += d*d; }
                *Y = 1.0 - dsm*den;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            const size_t K1 = (L==N1) ? 1 : K, K2 = (L==N2) ? 1 : K;
            const size_t I1 = (L==N1) ? 0 : B*(L-1), I2 = (L==N2) ? 0 : B*(L-1);
            for (size_t g=0; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    for (size_t l=0; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        XI1[l].val = *X1; XI1[l].ind = l;
                        XI2[l].val = *X2; XI2[l].ind = l;
                    }
                    qsort(XI1,L,sizeof(DBL),comp); qsort(XI2,L,sizeof(DBL),comp);
                    for (size_t l=0; l<L; ++l) { r1[XI1[l].ind] = r2[XI2[l].ind] = (double)l; }
                    dsm = 0.0; for (size_t l=0; l<L; ++l) { d = r1[l]-r2[l]; dsm += d*d; }
                    *Y = 1.0 - dsm*den;
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
