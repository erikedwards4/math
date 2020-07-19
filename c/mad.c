//Gets MAD (median absolute deviation from the median) along dim of X.

//The in-place versions still return the same Y, but modify X during processing.
//However, it turns out to be almost the identical speed for matrices.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mad_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mad_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int mad_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mad_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int mad_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mad_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t G = N / (B*L);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in mad_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (L==1)
    {
        const float z = 0.0f;
        cblas_scopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        cblas_scopy((int)L,X,1,X1,1);
        if (LAPACKE_slasrt_work('I',(int)N,X1)) { fprintf(stderr,"error in mad_s: problem with LAPACKE function\n"); }
        *Y = (N%2) ? X1[N/2] : 0.5f*(X1[N/2]+X1[N/2-1]);
        for (size_t l=0; l<L; l++) { X1[l] = fabsf(X1[l]-*Y); }
        if (LAPACKE_slasrt_work('I',(int)N,X1)) { fprintf(stderr,"error in mad_s: problem with LAPACKE function\n"); }
        *Y = (N%2) ? X1[N/2] : 0.5f*(X1[N/2]+X1[N/2-1]);
    }
    else
    {
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=J)
            {
                cblas_scopy((int)L,X,(int)K,X1,1);
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_s: problem with LAPACKE function\n"); }
                *Y++ = (L%2) ? X1[L/2] : 0.5f*(X1[L/2]+X1[L/2-1]);
                for (size_t l=0; l<L; l++) { X1[l] = fabsf(X1[l]-*Y); }
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_s: problem with LAPACKE function\n"); }
                *Y++ = (L%2) ? X1[L/2] : 0.5f*(X1[L/2]+X1[L/2-1]);
            }
        }
    }

    free(X1);
    return 0;
}


int mad_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mad_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t G = N / (B*L);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in mad_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (L==1)
    {
        const double z = 0.0;
        cblas_dcopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        cblas_dcopy((int)L,X,1,X1,1);
        if (LAPACKE_dlasrt_work('I',(int)N,X1)) { fprintf(stderr,"error in mad_d: problem with LAPACKE function\n"); }
        *Y = (N%2) ? X1[N/2] : 0.5*(X1[N/2]+X1[N/2-1]);
        for (size_t l=0; l<L; l++) { X1[l] = fabs(X1[l]-*Y); }
        if (LAPACKE_dlasrt_work('I',(int)N,X1)) { fprintf(stderr,"error in mad_d: problem with LAPACKE function\n"); }
        *Y = (N%2) ? X1[N/2] : 0.5*(X1[N/2]+X1[N/2-1]);
    }
    else
    {
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=J)
            {
                cblas_dcopy((int)L,X,(int)K,X1,1);
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_d: problem with LAPACKE function\n"); }
                *Y++ = (L%2) ? X1[L/2] : 0.5*(X1[L/2]+X1[L/2-1]);
                for (size_t l=0; l<L; l++) { X1[l] = fabs(X1[l]-*Y); }
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_d: problem with LAPACKE function\n"); }
                *Y++ = (L%2) ? X1[L/2] : 0.5*(X1[L/2]+X1[L/2-1]);
            }
        }
    }

    free(X1);
    return 0;
}


int mad_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mad_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t G = N / (B*L);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    if (L==1)
    {
        const float z = 0.0f;
        cblas_scopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        if (LAPACKE_slasrt_work('I',(int)N,X)) { fprintf(stderr,"error in mad_inplace_s: problem with LAPACKE function\n"); }
        *Y = (N%2) ? X[N/2] : 0.5f*(X[N/2]+X[N/2-1]);
        for (size_t n=0; n<N; n++) { X[n] = fabsf(X[n]-*Y); }
        if (LAPACKE_slasrt_work('I',(int)N,X)) { fprintf(stderr,"error in mad_inplace_s: problem with LAPACKE function\n"); }
        *Y = (N%2) ? X[N/2] : 0.5f*(X[N/2]+X[N/2-1]);
    }
    else if (K==1)
    {
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=J)
            {
                if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in mad_inplace_s: problem with LAPACKE function\n"); }
                *Y = (L%2) ? X[L/2] : 0.5f*(X[L/2]+X[L/2-1]);
                for (size_t l=0; l<L; l++) { X[l] = fabsf(X[l]-*Y); }
                if (LAPACKE_slasrt_work('I',(int)L,X)) { fprintf(stderr,"error in mad_inplace_s: problem with LAPACKE function\n"); }
                *Y++ = (L%2) ? X[L/2] : 0.5f*(X[L/2]+X[L/2-1]);
            }
        }
    }
    else
    {
        float *X1;
        if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in mad_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=J)
            {
                cblas_scopy((int)L,X,(int)K,X1,1);
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_inplace_s: problem with LAPACKE function\n"); }
                *Y = (L%2) ? X1[L/2] : 0.5f*(X1[L/2]+X1[L/2-1]);
                for (size_t l=0; l<L; l++) { X1[l] = fabsf(X1[l]-*Y); }
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_inplace_s: problem with LAPACKE function\n"); }
                *Y++ = (L%2) ? X1[L/2] : 0.5f*(X1[L/2]+X1[L/2-1]);
            }
        }
        free(X1);
    }

    return 0;
}


int mad_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mad_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t G = N / (B*L);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    if (L==1)
    {
        const double z = 0.0;
        cblas_dcopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)N,X)) { fprintf(stderr,"error in mad_inplace_d: problem with LAPACKE function\n"); }
        *Y = (N%2) ? X[N/2] : 0.5*(X[N/2]+X[N/2-1]);
        for (size_t n=0; n<N; n++) { X[n] = fabs(X[n]-*Y); }
        if (LAPACKE_dlasrt_work('I',(int)N,X)) { fprintf(stderr,"error in mad_inplace_d: problem with LAPACKE function\n"); }
        *Y = (N%2) ? X[N/2] : 0.5*(X[N/2]+X[N/2-1]);
    }
    else if (K==1)
    {
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=J)
            {
                if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in mad_inplace_d: problem with LAPACKE function\n"); }
                *Y = (L%2) ? X[L/2] : 0.5*(X[L/2]+X[L/2-1]);
                for (size_t l=0; l<L; l++) { X[l] = fabs(X[l]-*Y); }
                if (LAPACKE_dlasrt_work('I',(int)L,X)) { fprintf(stderr,"error in mad_inplace_d: problem with LAPACKE function\n"); }
                *Y++ = (L%2) ? X[L/2] : 0.5*(X[L/2]+X[L/2-1]);
            }
        }
    }
    else
    {
        double *X1;
        if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in mad_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=J)
            {
                cblas_dcopy((int)L,X,(int)K,X1,1);
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_inplace_d: problem with LAPACKE function\n"); }
                *Y = (L%2) ? X1[L/2] : 0.5*(X1[L/2]+X1[L/2-1]);
                for (size_t l=0; l<L; l++) { X1[l] = fabs(X1[l]-*Y); }
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mad_inplace_d: problem with LAPACKE function\n"); }
                *Y++ = (L%2) ? X1[L/2] : 0.5*(X1[L/2]+X1[L/2-1]);
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
