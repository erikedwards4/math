//Gets the coefficient of variation (std/mean) of each row or col of X according to dim.
//For complex case, output is real.
//This works in place.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int coeff_var_s (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased);
int coeff_var_d (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased);


int coeff_var_s (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t N2 = N/N1;
    const float z = 0.0f, o = 1.0f, ni = 1.0f / N1;
    const float den = (biased) ? 1.0f/N1 : 1.0f/(N1-1);

    if (N1==1) { cblas_scopy((int)N,&z,0,Y,1); }
    else if (N1==N)
    {
        Y[0] = 0.0f;
        for (size_t n=0; n<N; n++) { Y[0] += X[n]; }
        cblas_saxpy((int)N,-Y[0],&ni,0,X,1);
        Y[0] = sqrtf(cblas_sdot((int)N,X,1,X,1)*den) / Y[0];
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : (int)N1;
        const int ldb = (iscolmajor) ? (int)N2 : 1;
        const int ldc = (iscolmajor) ? (int)R : (int)C;
        float *x1;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)N1,&ni,0,x1,1);
        cblas_sgemv(Ord,Tr,(int)R,(int)C,1.0f,X,ldc,x1,1,0.0f,Y,1);
        cblas_scopy((int)N1,&o,0,x1,1);
        if (dim==0) { cblas_sgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,x1,lda,Y,ldb,1.0f,X,ldc); }
        else { cblas_sgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,Y,ldb,x1,lda,1.0f,X,ldc); }
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++) { Y[n2] = sqrtf(cblas_sdot((int)N1,&X[n2*N1],1,&X[n2*N1],1)*den) / Y[n2]; }
        }
        else
        {
            float *mn;
            if (!(mn=(float *)malloc((size_t)N2*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
            cblas_scopy((int)N2,Y,1,mn,1);
            for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
            cblas_sgemv(Ord,Tr,(int)R,(int)C,den,X,ldc,x1,1,0.0f,Y,1);
            for (size_t n2=0; n2<N2; n2++) { Y[n2] = sqrtf(Y[n2]) / mn[n2]; }
        }
        free(x1);
    }
    else if (iscolmajor && dim==0)
    {
        float *x1, *xni;
        if (!(x1=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)R,&o,0,x1,1); cblas_scopy((int)R,&ni,0,xni,1);
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, n+=RC, n2+=C)
            {
                cblas_sgemv(CblasColMajor,CblasTrans,(int)R,(int)C,1.0f,&X[n],(int)R,xni,1,0.0f,&Y[n2],1);
                cblas_sgemm(CblasColMajor,CblasTrans,CblasTrans,(int)R,(int)C,1,-1.0f,x1,1,&Y[n2],(int)C,1.0f,&X[n],(int)R);
                for (size_t c=0; c<C; c++) { Y[n2+c] = sqrtf(cblas_sdot((int)R,&X[n+c*R],1,&X[n+c*R],1)*den) / Y[n2+c]; }
            }
        }
        free(x1); free(xni);
    }
    else if (iscolmajor && dim==1)
    {
        float *x1, *xni;
        if (!(x1=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)C,&o,0,x1,1); cblas_scopy((int)C,&ni,0,xni,1);
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, n+=RC, n2+=R)
            {
                cblas_sgemv(CblasColMajor,CblasNoTrans,(int)R,(int)C,1.0f,&X[n],(int)R,xni,1,0.0f,&Y[n2],1);
                cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,1,-1.0f,&Y[n2],(int)R,x1,1,1.0f,&X[n],(int)R);
                for (size_t r=0; r<R; r++) { Y[n2+r] = sqrtf(cblas_sdot((int)C,&X[n+r],(int)R,&X[n+r],(int)R)*den) / Y[n2+r]; }
            }
        }
        free(x1); free(xni);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        float *x1, *xni;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)N1,&o,0,x1,1); cblas_scopy((int)N1,&ni,0,xni,1);
        for (size_t r=0, n=0, n2=0; r<R; r++, n+=C*SH, n2+=C)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0f,&X[n],(int)S,xni,1,0.0f,&Y[n2],1);
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)C,(int)S,1,-1.0f,&Y[n2],1,x1,(int)S,1.0f,&X[n],(int)S);
            for (size_t c=0; c<C; c++) { Y[n2+c] = sqrtf(cblas_sdot((int)N1,&X[n+c*S],1,&X[n+c*S],1)*den) / Y[n2+c]; }
        }
        free(x1); free(xni);
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, n+=J, n2++)
            {
                Y[n2] = cblas_sdot((int)N1,&X[n],(int)K,&ni,0);
                cblas_saxpy((int)N1,-Y[n2],&o,0,&X[n],(int)K);
                Y[n2] = sqrtf(cblas_sdot((int)N1,&X[n],(int)K,&X[n],(int)K)*den) / Y[n2];
            }
        }
    }
    
    return 0;
}


int coeff_var_d (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t N2 = N/N1;
    const double z = 0.0, o = 1.0, ni = 1.0 / N1;
    const double den = (biased) ? 1.0/N1 : 1.0/(N1-1);

    if (N1==1) { cblas_dcopy((int)N,&z,0,Y,1); }
    else if (N1==N)
    {
        Y[0] = 0.0;
        for (size_t n=0; n<N; n++) { Y[0] += X[n]; }
        cblas_daxpy((int)N,-Y[0],&ni,0,X,1);
        Y[0] = sqrt(cblas_ddot((int)N,X,1,X,1)*den) / Y[0];
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : (int)N1;
        const int ldb = (iscolmajor) ? (int)N2 : 1;
        const int ldc = (iscolmajor) ? (int)R : (int)C;
        double *x1;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)N1,&ni,0,x1,1);
        cblas_dgemv(Ord,Tr,(int)R,(int)C,1.0,X,ldc,x1,1,0.0,Y,1);
        cblas_dcopy((int)N1,&o,0,x1,1);
        if (dim==0) { cblas_dgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0,x1,lda,Y,ldb,1.0,X,ldc); }
        else { cblas_dgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0,Y,ldb,x1,lda,1.0,X,ldc); }
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n2=0; n2<N2; n2++) { Y[n2] = sqrt(cblas_ddot((int)N1,&X[n2*N1],1,&X[n2*N1],1)*den) / Y[n2]; }
        }
        else
        {
            double *mn;
            if (!(mn=(double *)malloc((size_t)N2*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
            cblas_dcopy((int)N2,Y,1,mn,1);
            for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
            cblas_dgemv(Ord,Tr,(int)R,(int)C,den,X,ldc,x1,1,0.0,Y,1);
            for (size_t n2=0; n2<N2; n2++) { Y[n2] = sqrt(Y[n2]) / mn[n2]; }
        }
        free(x1);
    }
    else if (iscolmajor && dim==0)
    {
        double *x1, *xni;
        if (!(x1=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)R,&o,0,x1,1); cblas_dcopy((int)R,&ni,0,xni,1);
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, n+=RC, n2+=C)
            {
                cblas_dgemv(CblasColMajor,CblasTrans,(int)R,(int)C,1.0,&X[n],(int)R,xni,1,0.0,&Y[n2],1);
                cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,(int)R,(int)C,1,-1.0,x1,1,&Y[n2],(int)C,1.0,&X[n],(int)R);
                for (size_t c=0; c<C; c++) { Y[n2+c] = sqrt(cblas_ddot((int)R,&X[n+c*R],1,&X[n+c*R],1)*den) / Y[n2+c]; }
            }
        }
        free(x1); free(xni);
    }
    else if (iscolmajor && dim==1)
    {
        double *x1, *xni;
        if (!(x1=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)C,&o,0,x1,1); cblas_dcopy((int)C,&ni,0,xni,1);
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, n+=RC, n2+=R)
            {
                cblas_dgemv(CblasColMajor,CblasNoTrans,(int)R,(int)C,1.0,&X[n],(int)R,xni,1,0.0,&Y[n2],1);
                cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,1,-1.0,&Y[n2],(int)R,x1,1,1.0,&X[n],(int)R);
                for (size_t r=0; r<R; r++) { Y[n2+r] = sqrt(cblas_ddot((int)C,&X[n+r],(int)R,&X[n+r],(int)R)*den) / Y[n2+r]; }
            }
        }
        free(x1); free(xni);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        double *x1, *xni;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)N1,&o,0,x1,1); cblas_dcopy((int)N1,&ni,0,xni,1);
        for (size_t r=0, n=0, n2=0; r<R; r++, n+=C*SH, n2+=C)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0,&X[n],(int)S,xni,1,0.0,&Y[n2],1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)C,(int)S,1,-1.0,&Y[n2],1,x1,(int)S,1.0,&X[n],(int)S);
            for (size_t c=0; c<C; c++) { Y[n2+c] = sqrt(cblas_ddot((int)N1,&X[n+c*S],1,&X[n+c*S],1)*den) / Y[n2+c]; }
        }
        free(x1); free(xni);
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (size_t l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, n+=J, n2++)
            {
                Y[n2] = cblas_ddot((int)N1,&X[n],(int)K,&ni,0);
                cblas_daxpy((int)N1,-Y[n2],&o,0,&X[n],(int)K);
                Y[n2] = sqrt(cblas_ddot((int)N1,&X[n],(int)K,&X[n],(int)K)*den) / Y[n2];
            }
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
