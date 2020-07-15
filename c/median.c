//Gets medians along dim of X.

//The in-place versions still return the same Y, but modify X during processing.
//However, it turns out to be almost the identical speed for matrices.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int median_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int median_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int median_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int median_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int median_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in median_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    float *X1;
    if (!(X1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in median_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N1==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        cblas_scopy((int)N1,X,1,X1,1);
        if (LAPACKE_slasrt_work('I',(int)N,X1)) { fprintf(stderr,"error in median_s: problem with LAPACKE function\n"); }
        *Y = (N%2) ? X1[N/2] : 0.5f*(X1[N/2]+X1[N/2-1]);
    }
    else
    {
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                cblas_scopy((int)N1,X,(int)K,X1,1);
                if (LAPACKE_slasrt_work('I',(int)N1,X1)) { fprintf(stderr,"error in median_s: problem with LAPACKE function\n"); }
                *Y++ = (N1%2) ? X1[N1/2] : 0.5f*(X1[N1/2]+X1[N1/2-1]);
            }
        }
    }

    free(X1);
    return 0;
}


int median_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in median_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    double *X1;
    if (!(X1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in median_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N1==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        cblas_dcopy((int)N1,X,1,X1,1);
        if (LAPACKE_dlasrt_work('I',(int)N,X1)) { fprintf(stderr,"error in median_d: problem with LAPACKE function\n"); }
        *Y = (N%2) ? X1[N/2] : 0.5*(X1[N/2]+X1[N/2-1]);
    }
    else
    {
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                cblas_dcopy((int)N1,X,(int)K,X1,1);
                if (LAPACKE_dlasrt_work('I',(int)N1,X1)) { fprintf(stderr,"error in median_d: problem with LAPACKE function\n"); }
                *Y++ = (N1%2) ? X1[N1/2] : 0.5*(X1[N1/2]+X1[N1/2-1]);
            }
        }
    }

    free(X1);
    return 0;
}


int median_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in median_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    if (N1==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        if (LAPACKE_slasrt_work('I',(int)N,X)) { fprintf(stderr,"error in median_inplace_s: problem with LAPACKE function\n"); }
        *Y = (N%2) ? X[N/2] : 0.5f*(X[N/2]+X[N/2-1]);
    }
    else if (K==1)
    {
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                if (LAPACKE_slasrt_work('I',(int)N1,X)) { fprintf(stderr,"error in median_inplace_s: problem with LAPACKE function\n"); }
                *Y++ = (N1%2) ? X[N1/2] : 0.5f*(X[N1/2]+X[N1/2-1]);
            }
        }
    }
    else
    {
        float *X1;
        if (!(X1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in median_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                cblas_scopy((int)N1,X,(int)K,X1,1);
                if (LAPACKE_slasrt_work('I',(int)N1,X1)) { fprintf(stderr,"error in median_inplace_s: problem with LAPACKE function\n"); }
                *Y++ = (N1%2) ? X1[N1/2] : 0.5f*(X1[N1/2]+X1[N1/2-1]);
            }
        }
        free(X1);
    }

    return 0;
}


int median_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in median_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    if (N1==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)N,X)) { fprintf(stderr,"error in median_inplace_d: problem with LAPACKE function\n"); }
        *Y = (N%2) ? X[N/2] : 0.5*(X[N/2]+X[N/2-1]);
    }
    else if (K==1)
    {
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                if (LAPACKE_dlasrt_work('I',(int)N1,X)) { fprintf(stderr,"error in median_inplace_d: problem with LAPACKE function\n"); }
                *Y++ = (N1%2) ? X[N1/2] : 0.5*(X[N1/2]+X[N1/2-1]);
            }
        }
    }
    else
    {
        double *X1;
        if (!(X1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in median_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                cblas_dcopy((int)N1,X,(int)K,X1,1);
                if (LAPACKE_dlasrt_work('I',(int)N1,X1)) { fprintf(stderr,"error in median_inplace_d: problem with LAPACKE function\n"); }
                *Y++ = (N1%2) ? X1[N1/2] : 0.5*(X1[N1/2]+X1[N1/2-1]);
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
