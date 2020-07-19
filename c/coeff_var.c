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

int coeff_var_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int coeff_var_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);


int coeff_var_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in coeff_var_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t V = N/L;
    const float z = 0.0f, o = 1.0f, ni = 1.0f / L;
    const float den = (biased) ? 1.0f/L : 1.0f/(L-1);

    if (L==1) { cblas_scopy((int)N,&z,0,Y,1); }
    else if (L==N)
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
        const int lda = (iscolmajor) ? 1 : (int)L;
        const int ldb = (iscolmajor) ? (int)V : 1;
        const int ldc = (iscolmajor) ? (int)R : (int)C;
        float *x1;
        if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)L,&ni,0,x1,1);
        cblas_sgemv(Ord,Tr,(int)R,(int)C,1.0f,X,ldc,x1,1,0.0f,Y,1);
        cblas_scopy((int)L,&o,0,x1,1);
        if (dim==0) { cblas_sgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,x1,lda,Y,ldb,1.0f,X,ldc); }
        else { cblas_sgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,Y,ldb,x1,lda,1.0f,X,ldc); }
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t v=0; v<V; v++, X+=L) { Y[v] = sqrtf(cblas_sdot((int)L,X,1,X,1)*den) / Y[v]; }
        }
        else
        {
            float *mn;
            if (!(mn=(float *)malloc(V*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
            cblas_scopy((int)V,Y,1,mn,1);
            for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
            cblas_sgemv(Ord,Tr,(int)R,(int)C,den,X,ldc,x1,1,0.0f,Y,1);
            for (size_t v=0; v<V; v++) { Y[v] = sqrtf(Y[v]) / mn[v]; }
        }
        free(x1);
    }
    else if (iscolmajor && dim==0)
    {
        float *x1, *xni;
        if (!(x1=(float *)malloc(R*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc(R*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)R,&o,0,x1,1); cblas_scopy((int)R,&ni,0,xni,1);
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, X+=RC, Y+=C)
            {
                cblas_sgemv(CblasColMajor,CblasTrans,(int)R,(int)C,1.0f,X,(int)R,xni,1,0.0f,Y,1);
                cblas_sgemm(CblasColMajor,CblasTrans,CblasTrans,(int)R,(int)C,1,-1.0f,x1,1,Y,(int)C,1.0f,X,(int)R);
                for (size_t c=0; c<C; c++, X+=R, Y++) { *Y = sqrtf(cblas_sdot((int)R,X,1,X,1)*den) / *Y; }
            }
        }
        free(x1); free(xni);
    }
    else if (iscolmajor && dim==1)
    {
        float *x1, *xni;
        if (!(x1=(float *)malloc(C*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc(C*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)C,&o,0,x1,1); cblas_scopy((int)C,&ni,0,xni,1);
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, X+=(R-1)*C)
            {
                cblas_sgemv(CblasColMajor,CblasNoTrans,(int)R,(int)C,1.0f,X,(int)R,xni,1,0.0f,Y,1);
                cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,1,-1.0f,Y,(int)R,x1,1,1.0f,X,(int)R);
                for (size_t r=0; r<R; r++, X++, Y++) { *Y = sqrtf(cblas_sdot((int)C,X,(int)R,X,(int)R)*den) / *Y; }
            }
        }
        free(x1); free(xni);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        float *x1, *xni;
        if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)L,&o,0,x1,1); cblas_scopy((int)L,&ni,0,xni,1);
        for (size_t r=0; r<R; r++)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0f,X,(int)S,xni,1,0.0f,Y,1);
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)C,(int)S,1,-1.0f,Y,1,x1,(int)S,1.0f,X,(int)S);
            for (size_t c=0; c<C; c++, X+=S, Y++) { *Y = sqrtf(cblas_sdot((int)L,X,1,X,1)*den) / *Y; }
        }
        free(x1); free(xni);
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
                cblas_saxpy((int)L,-*Y,&o,0,X,(int)K);
                *Y = sqrtf(cblas_sdot((int)L,X,(int)K,X,(int)K)*den) / *Y;
            }
        }
    }
    
    return 0;
}


int coeff_var_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in coeff_var_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t V = N/L;
    const double z = 0.0, o = 1.0, ni = 1.0 / L;
    const double den = (biased) ? 1.0/L : 1.0/(L-1);

    if (L==1) { cblas_dcopy((int)N,&z,0,Y,1); }
    else if (L==N)
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
        const int lda = (iscolmajor) ? 1 : (int)L;
        const int ldb = (iscolmajor) ? (int)V : 1;
        const int ldc = (iscolmajor) ? (int)R : (int)C;
        double *x1;
        if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)L,&ni,0,x1,1);
        cblas_dgemv(Ord,Tr,(int)R,(int)C,1.0,X,ldc,x1,1,0.0,Y,1);
        cblas_dcopy((int)L,&o,0,x1,1);
        if (dim==0) { cblas_dgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0,x1,lda,Y,ldb,1.0,X,ldc); }
        else { cblas_dgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0,Y,ldb,x1,lda,1.0,X,ldc); }
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t v=0; v<V; v++, X+=L) { Y[v] = sqrt(cblas_ddot((int)L,X,1,X,1)*den) / Y[v]; }
        }
        else
        {
            double *mn;
            if (!(mn=(double *)malloc(V*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
            cblas_dcopy((int)V,Y,1,mn,1);
            for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
            cblas_dgemv(Ord,Tr,(int)R,(int)C,den,X,ldc,x1,1,0.0,Y,1);
            for (size_t v=0; v<V; v++) { Y[v] = sqrt(Y[v]) / mn[v]; }
        }
        free(x1);
    }
    else if (iscolmajor && dim==0)
    {
        double *x1, *xni;
        if (!(x1=(double *)malloc(R*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc(R*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)R,&o,0,x1,1); cblas_dcopy((int)R,&ni,0,xni,1);
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, X+=RC, Y+=C)
            {
                cblas_dgemv(CblasColMajor,CblasTrans,(int)R,(int)C,1.0,X,(int)R,xni,1,0.0,Y,1);
                cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,(int)R,(int)C,1,-1.0,x1,1,Y,(int)C,1.0,X,(int)R);
                for (size_t c=0; c<C; c++, X+=R, Y++) { *Y = sqrt(cblas_ddot((int)R,X,1,X,1)*den) / *Y; }
            }
        }
        free(x1); free(xni);
    }
    else if (iscolmajor && dim==1)
    {
        double *x1, *xni;
        if (!(x1=(double *)malloc(C*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc(C*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)C,&o,0,x1,1); cblas_dcopy((int)C,&ni,0,xni,1);
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, X+=(R-1)*C)
            {
                cblas_dgemv(CblasColMajor,CblasNoTrans,(int)R,(int)C,1.0,X,(int)R,xni,1,0.0,Y,1);
                cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,1,-1.0,Y,(int)R,x1,1,1.0,X,(int)R);
                for (size_t r=0; r<R; r++, X++, Y++) { *Y = sqrt(cblas_ddot((int)C,X,(int)R,X,(int)R)*den) / *Y; }
            }
        }
        free(x1); free(xni);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        double *x1, *xni;
        if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)L,&o,0,x1,1); cblas_dcopy((int)L,&ni,0,xni,1);
        for (size_t r=0; r<R; r++)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0,X,(int)S,xni,1,0.0,Y,1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)C,(int)S,1,-1.0,Y,1,x1,(int)S,1.0,X,(int)S);
            for (size_t c=0; c<C; c++, X+=S, Y++) { *Y = sqrt(cblas_ddot((int)L,X,1,X,1)*den) / *Y; }
        }
        free(x1); free(xni);
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
                cblas_daxpy((int)L,-*Y,&o,0,X,(int)K);
                *Y = sqrt(cblas_ddot((int)L,X,(int)K,X,(int)K)*den) / *Y;
            }
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
