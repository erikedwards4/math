//Gets all percentiles in P along dim of X.
//P is of length Q, so Y ends up with length Q along dim.

//The in-place versions still return the same Y, but modify X during processing.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int prctiles_s (float *Y, const float *X, const float *P, const size_t Q, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int prctiles_d (double *Y, const double *X, const double *P, const size_t Q, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);

int prctiles_inplace_s (float *Y, float *X, const float *P, const size_t Q, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int prctiles_inplace_d (double *Y, double *X, const double *P, const size_t Q, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);


int prctiles_s (float *Y, const float *X, const float *P, const size_t Q, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    float *X1, *w1, *w2;
    size_t *i1, *i2;
    if (!(X1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in prctiles_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w1=(float *)malloc(Q*sizeof(float)))) { fprintf(stderr,"error in prctiles_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w2=(float *)malloc(Q*sizeof(float)))) { fprintf(stderr,"error in prctiles_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(i1=(size_t *)malloc(Q*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(i2=(size_t *)malloc(Q*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_s: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t q=0; q<Q; q++)
    {
        if (P[q]<0.0f || P[q]>100.0f) { fprintf(stderr,"error in prctiles_s: prctiles must be in [0 100]"); return 1; }
        float p1 = (P[q]/100.0f)*(N1-1);
        i1[q] = (P[q]<100.0f) ? (size_t)floorf(p1) : N1-2;
        i2[q] = i1[q] + 1;
        w2[q] = (P[q]<100.0f) ? p1-floorf(p1) : 1.0f;
        w1[q] = 1.0f - w2[q];
    }
    
    if (N1==1 && Q==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        cblas_scopy((int)N1,X,1,X1,1);
        if (LAPACKE_slasrt_work('I',(int)N,X1)) { fprintf(stderr,"error in prctiles_s: problem with LAPACKE function\n"); }
        for (size_t q=0; q<Q; q++) { *Y++ = w1[q]*X1[i1[q]] + w2[q]*X1[i2[q]]; }
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const size_t Jy = (iscolmajor) ? ((dim==0) ? Q : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : Q);
        for (size_t l=0; l<L; l++, X+=M*(N1-J), Y+=M*(Q-Jy))
        {
            for (size_t m=0; m<M; m++, X+=J, Y-=Q*K-Jy)
            {
                cblas_scopy((int)N1,X,(int)K,X1,1);
                if (LAPACKE_slasrt_work('I',(int)N1,X1)) { fprintf(stderr,"error in prctiles_s: problem with LAPACKE function\n"); }
                for (size_t q=0; q<Q; q++, Y+=K) { *Y = w1[q]*X1[i1[q]] + w2[q]*X1[i2[q]]; }
            }
        }
    }

    free(X1);
    return 0;
}


int prctiles_d (double *Y, const double *X, const double *P, const size_t Q, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    double *X1, *w1, *w2;
    size_t *i1, *i2;
    if (!(X1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in prctiles_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w1=(double *)malloc(Q*sizeof(double)))) { fprintf(stderr,"error in prctiles_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w2=(double *)malloc(Q*sizeof(double)))) { fprintf(stderr,"error in prctiles_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(i1=(size_t *)malloc(Q*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(i2=(size_t *)malloc(Q*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_d: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t q=0; q<Q; q++)
    {
        if (P[q]<0.0 || P[q]>100.0) { fprintf(stderr,"error in prctiles_d: prctiles must be in [0 100]"); return 1; }
        double p1 = (P[q]/100.0)*(N1-1);
        i1[q] = (P[q]<100.0) ? (size_t)floor(p1) : N1-2;
        i2[q] = i1[q] + 1;
        w2[q] = (P[q]<100.0) ? p1-floor(p1) : 1.0;
        w1[q] = 1.0 - w2[q];
    }
    
    if (N1==1 && Q==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        cblas_dcopy((int)N1,X,1,X1,1);
        if (LAPACKE_dlasrt_work('I',(int)N,X1)) { fprintf(stderr,"error in prctiles_d: problem with LAPACKE function\n"); }
        for (size_t q=0; q<Q; q++) { *Y++ = w1[q]*X1[i1[q]] + w2[q]*X1[i2[q]]; }
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const size_t Jy = (iscolmajor) ? ((dim==0) ? Q : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : Q);
        for (size_t l=0; l<L; l++, X+=M*(N1-J), Y+=M*(Q-Jy))
        {
            for (size_t m=0; m<M; m++, X+=J, Y-=Q*K-Jy)
            {
                cblas_dcopy((int)N1,X,(int)K,X1,1);
                if (LAPACKE_dlasrt_work('I',(int)N1,X1)) { fprintf(stderr,"error in prctiles_d: problem with LAPACKE function\n"); }
                for (size_t q=0; q<Q; q++, Y+=K) { *Y = w1[q]*X1[i1[q]] + w2[q]*X1[i2[q]]; }
            }
        }
    }

    free(X1);
    return 0;
}


int prctiles_inplace_s (float *Y, float *X, const float *P, const size_t Q, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    float *w1, *w2;
    size_t *i1, *i2;
    if (!(w1=(float *)malloc(Q*sizeof(float)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w2=(float *)malloc(Q*sizeof(float)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(i1=(size_t *)malloc(Q*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(i2=(size_t *)malloc(Q*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t q=0; q<Q; q++)
    {
        if (P[q]<0.0f || P[q]>100.0f) { fprintf(stderr,"error in prctiles_inplace_s: prctiles must be in [0 100]"); return 1; }
        float p1 = (P[q]/100.0f)*(N1-1);
        i1[q] = (P[q]<100.0f) ? (size_t)floorf(p1) : N1-2;
        i2[q] = i1[q] + 1;
        w2[q] = (P[q]<100.0f) ? p1-floorf(p1) : 1.0f;
        w1[q] = 1.0f - w2[q];
    }
    
    if (N1==1 && Q==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        if (LAPACKE_slasrt_work('I',(int)N,X)) { fprintf(stderr,"error in prctiles_inplace_s: problem with LAPACKE function\n"); }
        for (size_t q=0; q<Q; q++) { *Y++ = w1[q]*X[i1[q]] + w2[q]*X[i2[q]]; }
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const size_t Jy = (iscolmajor) ? ((dim==0) ? Q : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : Q);
        float *X1;
        if (!(X1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in prctiles_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++, X+=M*(N1-J), Y+=M*(Q-Jy))
        {
            for (size_t m=0; m<M; m++, X+=J, Y-=Q*K-Jy)
            {
                cblas_scopy((int)N1,X,(int)K,X1,1);
                if (LAPACKE_slasrt_work('I',(int)N1,X1)) { fprintf(stderr,"error in prctiles_inplace_s: problem with LAPACKE function\n"); }
                for (size_t q=0; q<Q; q++, Y+=K) { *Y = w1[q]*X1[i1[q]] + w2[q]*X1[i2[q]]; }
            }
        }
        free(X1);
    }

    return 0;
}


int prctiles_inplace_d (double *Y, double *X, const double *P, const size_t Q, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Prep interpolation
    double *w1, *w2;
    size_t *i1, *i2;
    if (!(w1=(double *)malloc(Q*sizeof(double)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(w2=(double *)malloc(Q*sizeof(double)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(i1=(size_t *)malloc(Q*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(i2=(size_t *)malloc(Q*sizeof(size_t)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
    for (size_t q=0; q<Q; q++)
    {
        if (P[q]<0.0 || P[q]>100.0) { fprintf(stderr,"error in prctiles_inplace_d: prctiles must be in [0 100]"); return 1; }
        double p1 = (P[q]/100.0)*(N1-1);
        i1[q] = (P[q]<100.0) ? (size_t)floor(p1) : N1-2;
        i2[q] = i1[q] + 1;
        w2[q] = (P[q]<100.0) ? p1-floor(p1) : 1.0;
        w1[q] = 1.0 - w2[q];
    }
    
    if (N1==1 && Q==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)N,X)) { fprintf(stderr,"error in prctiles_inplace_d: problem with LAPACKE function\n"); }
        for (size_t q=0; q<Q; q++) { *Y++ = w1[q]*X[i1[q]] + w2[q]*X[i2[q]]; }
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const size_t Jy = (iscolmajor) ? ((dim==0) ? Q : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : Q);
        double *X1;
        if (!(X1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in prctiles_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++, X+=M*(N1-J), Y+=M*(Q-Jy))
        {
            for (size_t m=0; m<M; m++, X+=J, Y-=Q*K-Jy)
            {
                cblas_dcopy((int)N1,X,(int)K,X1,1);
                if (LAPACKE_dlasrt_work('I',(int)N1,X1)) { fprintf(stderr,"error in prctiles_inplace_d: problem with LAPACKE function\n"); }
                for (size_t q=0; q<Q; q++, Y+=K) { *Y = w1[q]*X1[i1[q]] + w2[q]*X1[i2[q]]; }
            }
        }
        free(X1);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
