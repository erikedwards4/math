//Gets mean of each row or col of X according to dim.
//For complex case, real and imag parts calculated separately.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int mean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ni = 1.0f / L;

    if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        Y[0] = 0.0f;
        for (size_t n=0; n<N; n++) { Y[0] += X[n]; }
        Y[0] *= ni;
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        float *xni;
        if (!(xni=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in mean_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)L,&ni,0,xni,1);
        cblas_sgemv(Ord,Tr,(int)R,(int)C,1.0f,X,lda,xni,1,0.0f,Y,1);
        free(xni);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        float *xni;
        if (!(xni=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in mean_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)L,&ni,0,xni,1);
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, X+=RC, Y+=RC/L)
            {
                cblas_sgemv(CblasColMajor,Tr,(int)R,(int)C,1.0f,X,(int)R,xni,1,0.0f,Y,1);
            }
        }
        free(xni);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        float *xni;
        if (!(xni=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in mean_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)L,&ni,0,xni,1);
        for (size_t r=0; r<R; r++, X+=C*SH, Y+=C)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0f,X,(int)S,xni,1,0.0f,Y,1);
        }
        free(xni);
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=J, Y++)
            {
                *Y = cblas_sdot((int)L,X,(int)K,&ni,0);
            }
        }
    }
    
    return 0;
}


int mean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni = 1.0 / L;

    if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        Y[0] = 0.0;
        for (size_t n=0; n<N; n++) { Y[0] += X[n]; }
        Y[0] *= ni;
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        double *xni;
        if (!(xni=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in mean_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)L,&ni,0,xni,1);
        cblas_dgemv(Ord,Tr,(int)R,(int)C,1.0,X,lda,xni,1,0.0,Y,1);
        free(xni);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        double *xni;
        if (!(xni=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in mean_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)L,&ni,0,xni,1);
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, X+=RC, Y+=RC/L)
            {
                cblas_dgemv(CblasColMajor,Tr,(int)R,(int)C,1.0,X,(int)R,xni,1,0.0,Y,1);
            }
        }
        free(xni);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        double *xni;
        if (!(xni=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in mean_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)L,&ni,0,xni,1);
        for (size_t r=0; r<R; r++, X+=C*SH, Y+=C)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0,X,(int)S,xni,1,0.0,Y,1);
        }
        free(xni);
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=J, Y++)
            {
                *Y = cblas_ddot((int)L,X,(int)K,&ni,0);
            }
        }
    }
    
    return 0;
}


int mean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean_c: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ni[2] = {1.0f/L,0.0f};

    if (L==1) { cblas_ccopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        Y[0] = Y[1] = 0.0f;
        for (size_t n=0; n<2*L; n+=2) { Y[0] += X[n]; Y[1] += X[n+1]; }
        Y[0] *= ni[0]; Y[1] *= ni[0];
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
        float *xni;
        if (!(xni=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in mean_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy((int)L,ni,0,xni,1);
        cblas_cgemv(Ord,Tr,(int)R,(int)C,o,X,lda,xni,1,z,Y,1);
        free(xni);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
        float *xni;
        if (!(xni=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in mean_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy((int)L,o,0,xni,1);
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, X+=2*RC, Y+=2*RC/L)
            {
                cblas_cgemv(CblasColMajor,Tr,(int)R,(int)C,o,X,(int)R,xni,1,z,Y,1);
            }
        }
        free(xni);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
        float *xni;
        if (!(xni=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in mean_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy((int)L,o,0,xni,1);
        for (size_t r=0; r<R; r++, X+=2*C*SH, Y+=2*C)
        {
            cblas_cgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,o,X,(int)S,xni,1,z,Y,1);
        }
        free(xni);
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t g=0; g<G; g++, X+=2*B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=2*J, Y+=2)
            {
                cblas_cdotu_sub((int)L,X,(int)K,ni,0,(_Complex float *)Y);
            }
        }
    }
    
    return 0;
}


int mean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean_z: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni[2] = {1.0/L,0.0};

    if (L==1) { cblas_zcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        Y[0] = Y[1] = 0.0;
        for (size_t n=0; n<2*L; n+=2) { Y[0] += X[n]; Y[1] += X[n+1]; }
        Y[0] *= ni[0]; Y[1] *= ni[0];
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
        double *xni;
        if (!(xni=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in mean_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy((int)L,ni,0,xni,1);
        cblas_zgemv(Ord,Tr,(int)R,(int)C,o,X,lda,xni,1,z,Y,1);
        free(xni);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
        double *xni;
        if (!(xni=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in mean_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy((int)L,o,0,xni,1);
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, X+=2*RC, Y+=2*RC/L)
            {
                cblas_zgemv(CblasColMajor,Tr,(int)R,(int)C,o,X,(int)R,xni,1,z,Y,1);
            }
        }
        free(xni);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
        double *xni;
        if (!(xni=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in mean_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy((int)L,o,0,xni,1);
        for (size_t r=0; r<R; r++, X+=2*C*SH, Y+=2*C)
        {
            cblas_zgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,o,X,(int)S,xni,1,z,Y,1);
        }
        free(xni);
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t g=0; g<G; g++, X+=2*B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=2*J, Y+=2)
            {
                cblas_zdotu_sub((int)L,X,(int)K,ni,0,(_Complex double *)Y);
            }
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
