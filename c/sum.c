//Gets sum of each row or col of X according to dim.
//For complex case, real and imag parts calculated separately.

//For 2D case, this is definitely faster using cblas_?gemv.

//No in-place version, since cblas_?gemv doesn't work for that.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sum_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor);
int sum_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor);
int sum_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor);
int sum_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor);


int sum_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if (N1==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        Y[0] = 0.0f;
        for (size_t n=0; n<N; n++) { Y[0] += X[n]; }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const float o = 1.0f;
        float *x1;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in sum_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)N1,&o,0,x1,1);
        cblas_sgemv(Ord,Tr,(int)R,(int)C,1.0f,X,lda,x1,1,0.0f,Y,1);
        //for (size_t c=0; c<C; c++) { Y[c] = cblas_sdot((int)R,&X[c*R],1,&o,0); }
        free(x1);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const float o = 1.0f;
        float *x1;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in sum_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)N1,&o,0,x1,1);
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, n+=RC, n2+=RC/N1)
            {
                cblas_sgemv(CblasColMajor,Tr,(int)R,(int)C,1.0f,&X[n],(int)R,x1,1,0.0f,&Y[n2],1);
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        const float o = 1.0f;
        float *x1;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in sum_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)N1,&o,0,x1,1);
        for (size_t r=0, n=0, n2=0; r<R; r++, n+=C*SH, n2+=C)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0f,&X[n],(int)S,x1,1,0.0f,&Y[n2],1);
        }
        free(x1);
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const float o = 1.0f;
        for (size_t l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, n+=J, n2++)
            {
                Y[n2] = cblas_sdot((int)N1,&X[n],(int)K,&o,0);
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int sum_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        Y[0] = 0.0;
        for (size_t n=0; n<N; n++) { Y[0] += X[n]; }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const double o = 1.0;
        double *x1;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in sum_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)N1,&o,0,x1,1);
        cblas_dgemv(Ord,Tr,(int)R,(int)C,1.0,X,lda,x1,1,0.0,Y,1);
        free(x1);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const double o = 1.0;
        double *x1;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in sum_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)N1,&o,0,x1,1);
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, n+=RC, n2+=RC/N1)
            {
                cblas_dgemv(CblasColMajor,Tr,(int)R,(int)C,1.0,&X[n],(int)R,x1,1,0.0,&Y[n2],1);
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        const double o = 1.0;
        double *x1;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in sum_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)N1,&o,0,x1,1);
        for (size_t r=0, n=0, n2=0; r<R; r++, n+=C*SH, n2+=C)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0,&X[n],(int)S,x1,1,0.0,&Y[n2],1);
        }
        free(x1);
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const double o = 1.0;
        for (size_t l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, n+=J, n2++)
            {
                Y[n2] = cblas_ddot((int)N1,&X[n],(int)K,&o,0);
            }
        }
    }
    
    return 0;
}


int sum_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_ccopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        Y[0] = Y[1] = 0.0f;
        for (size_t n=0; n<2*N1; n+=2) { Y[0] += X[n]; Y[1] += X[n+1]; }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
        float *x1;
        if (!(x1=(float *)malloc((size_t)(2*N1)*sizeof(float)))) { fprintf(stderr,"error in sum_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy((int)N1,o,0,x1,1);
        cblas_cgemv(Ord,Tr,(int)R,(int)C,o,X,lda,x1,1,z,Y,1);
        free(x1);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
        float *x1;
        if (!(x1=(float *)malloc((size_t)(2*N1)*sizeof(float)))) { fprintf(stderr,"error in sum_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy((int)N1,o,0,x1,1);
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, n+=2*RC, n2+=2*RC/N1)
            {
                cblas_cgemv(CblasColMajor,Tr,(int)R,(int)C,o,&X[n],(int)R,x1,1,z,&Y[n2],1);
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
        float *x1;
        if (!(x1=(float *)malloc((size_t)(2*N1)*sizeof(float)))) { fprintf(stderr,"error in sum_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy((int)N1,o,0,x1,1);
        for (size_t r=0, n=0, n2=0; r<R; r++, n+=2*C*SH, n2+=2*C)
        {
            cblas_cgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,o,&X[n],(int)S,x1,1,z,&Y[n2],1);
        }
        free(x1);
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const float o[2] = {1.0f,0.0f};
        for (size_t l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
        {
            for (size_t m=0; m<M; m++, n+=2*J, n2+=2)
            {
                cblas_cdotu_sub((int)N1,&X[n],(int)K,o,0,(_Complex float *)(&Y[n2]));
            }
        }
    }
    
    return 0;
}


int sum_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_zcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        Y[0] = Y[1] = 0.0;
        for (size_t n=0; n<2*N1; n+=2) { Y[0] += X[n]; Y[1] += X[n+1]; }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? (int)R : (int)C;
        const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
        double *x1;
        if (!(x1=(double *)malloc((size_t)(2*N1)*sizeof(double)))) { fprintf(stderr,"error in sum_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy((int)N1,o,0,x1,1);
        cblas_zgemv(Ord,Tr,(int)R,(int)C,o,X,lda,x1,1,z,Y,1);
        free(x1);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
        double *x1;
        if (!(x1=(double *)malloc((size_t)(2*N1)*sizeof(double)))) { fprintf(stderr,"error in sum_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy((int)N1,o,0,x1,1);
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, n+=2*RC, n2+=2*RC/N1)
            {
                cblas_zgemv(CblasColMajor,Tr,(int)R,(int)C,o,&X[n],(int)R,x1,1,z,&Y[n2],1);
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
        double *x1;
        if (!(x1=(double *)malloc((size_t)(2*N1)*sizeof(double)))) { fprintf(stderr,"error in sum_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy((int)N1,o,0,x1,1);
        for (size_t r=0, n=0, n2=0; r<R; r++, n+=2*C*SH, n2+=2*C)
        {
            cblas_zgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,o,&X[n],(int)S,x1,1,z,&Y[n2],1);
        }
        free(x1);
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const double o[2] = {1.0,0.0};
        for (size_t l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
        {
            for (size_t m=0; m<M; m++, n+=2*J, n2+=2)
            {
                cblas_zdotu_sub((int)N1,&X[n],(int)K,o,0,(_Complex double *)(&Y[n2]));
            }
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
