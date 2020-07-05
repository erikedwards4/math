//Gets variance of each row or col of X according to dim.
//For complex case, output is real.
//This works in place.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int var_s (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased);
int var_d (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased);
int var_c (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased);
int var_z (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased);


int var_s (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    if (R<1) { fprintf(stderr,"error in var_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in var_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in var_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in var_s: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int N2 = N/N1;
    const float z = 0.0f, o = 1.0f, ni = 1.0f / N1;
    const float den = (biased) ? 1.0f/N1 : 1.0f/(N1-1);

    if (N1==1) { cblas_scopy(N,&z,0,Y,1); }
    else if (N1==N)
    {
        Y[0] = 0.0f;
        for (int n=0; n<N; n++) { Y[0] += X[n]; }
        cblas_saxpy(N,-Y[0],&ni,0,X,1);
        Y[0] = cblas_sdot(N,X,1,X,1) * den;
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : N1;
        const int ldb = (iscolmajor) ? RC/N1 : 1;
        const int ldc = (iscolmajor) ? R : C;
        float *x1;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N1,&ni,0,x1,1);
        cblas_sgemv(Ord,Tr,(int)R,(int)C,1.0f,X,ldc,x1,1,0.0f,Y,1);
        cblas_scopy(N1,&o,0,x1,1);
        if (dim==0) { cblas_sgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,x1,lda,Y,ldb,1.0f,X,ldc); }
        else { cblas_sgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,Y,ldb,x1,lda,1.0f,X,ldc); }
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (int n2=0; n2<N2; n2++) { Y[n2] = cblas_sdot(N1,&X[n2*N1],1,&X[n2*N1],1) * den; }
        }
        else
        {
            for (int n=0; n<N; n++) { X[n] *= X[n]; }
            cblas_sgemv(Ord,Tr,(int)R,(int)C,den,X,ldc,x1,1,0.0f,Y,1);
        }
        free(x1);
    }
    else if (iscolmajor && dim==0)
    {
        float *x1, *xni;
        if (!(x1=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)R,&o,0,x1,1); cblas_scopy((int)R,&ni,0,xni,1);
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, n+=RC, n2+=C)
            {
                cblas_sgemv(CblasColMajor,CblasTrans,(int)R,(int)C,1.0f,&X[n],(int)R,xni,1,0.0f,&Y[n2],1);
                cblas_sgemm(CblasColMajor,CblasTrans,CblasTrans,(int)R,(int)C,1,-1.0f,x1,1,&Y[n2],(int)C,1.0f,&X[n],(int)R);
                for (size_t c=0; c<C; c++) { Y[n2+c] = cblas_sdot((int)R,&X[n+c*R],1,&X[n+c*R],1) * den; }
            }
        }
        free(x1); free(xni);
    }
    else if (iscolmajor && dim==1)
    {
        float *x1, *xni;
        if (!(x1=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)C,&o,0,x1,1); cblas_scopy((int)C,&ni,0,xni,1);
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, n+=RC, n2+=R)
            {
                cblas_sgemv(CblasColMajor,CblasNoTrans,(int)R,(int)C,1.0f,&X[n],(int)R,xni,1,0.0f,&Y[n2],1);
                cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,1,-1.0f,&Y[n2],(int)R,x1,1,1.0f,&X[n],(int)R);
                for (size_t r=0; r<R; r++) { Y[n2+r] = cblas_sdot((int)C,&X[n+r],(int)R,&X[n+r],(int)R) * den; }
            }
        }
        free(x1); free(xni);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        float *x1, *xni;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N1,&o,0,x1,1); cblas_scopy(N1,&ni,0,xni,1);
        for (size_t r=0, n=0, n2=0; r<R; r++, n+=C*SH, n2+=C)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0f,&X[n],(int)S,xni,1,0.0f,&Y[n2],1);
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)C,(int)S,1,-1.0f,&Y[n2],1,x1,(int)S,1.0f,&X[n],(int)S);
            for (size_t c=0; c<C; c++) { Y[n2+c] = cblas_sdot(N1,&X[n+c*S],1,&X[n+c*S],1) * den; }
        }
        free(x1); free(xni);
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (int l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
        {
            for (int m=0; m<M; m++, n+=J, n2++)
            {
                Y[n2] = cblas_sdot(N1,&X[n],K,&ni,0);
                cblas_saxpy(N1,-Y[n2],&o,0,&X[n],K);
                Y[n2] = cblas_sdot(N1,&X[n],K,&X[n],K) * den;
            }
        }
    }
    
    return 0;
}


int var_d (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    if (R<1) { fprintf(stderr,"error in var_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in var_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in var_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in var_d: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double z = 0.0, o = 1.0, ni = 1.0 / N1;
    const double den = (biased) ? 1.0/N1 : 1.0/(N1-1);

    if (N1==1) { cblas_dcopy(N,&z,0,Y,1); }
    else if (N1==N)
    {
        Y[0] = 0.0;
        for (int n=0; n<N; n++) { Y[0] += X[n]; }
        cblas_daxpy(N,-Y[0],&ni,0,X,1);
        Y[0] = cblas_ddot(N,X,1,X,1) * den;
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : N1;
        const int ldb = (iscolmajor) ? RC/N1 : 1;
        const int ldc = (iscolmajor) ? R : C;
        double *x1;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N1,&ni,0,x1,1);
        cblas_dgemv(Ord,Tr,(int)R,(int)C,1.0,X,ldc,x1,1,0.0,Y,1);
        cblas_dcopy(N1,&o,0,x1,1);
        if (dim==0) { cblas_dgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,x1,lda,Y,ldb,1.0f,X,ldc); }
        else { cblas_dgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,Y,ldb,x1,lda,1.0f,X,ldc); }
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (int n2=0; n2<N2; n2++) { Y[n2] = cblas_ddot(N1,&X[n2*N1],1,&X[n2*N1],1) * den; }
        }
        else
        {
            for (int n=0; n<N; n++) { X[n] *= X[n]; }
            cblas_dgemv(Ord,Tr,(int)R,(int)C,den,X,ldc,x1,1,0.0,Y,1);
        }
        free(x1);
    }
    else if (iscolmajor && dim==0)
    {
        double *x1, *xni;
        if (!(x1=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)R,&o,0,x1,1); cblas_dcopy((int)R,&ni,0,xni,1);
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, n+=RC, n2+=C)
            {
                cblas_dgemv(CblasColMajor,CblasTrans,(int)R,(int)C,1.0,&X[n],(int)R,xni,1,0.0,&Y[n2],1);
                cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,(int)R,(int)C,1,-1.0,x1,1,&Y[n2],(int)C,1.0,&X[n],(int)R);
                for (size_t c=0; c<C; c++) { Y[n2+c] = cblas_ddot((int)R,&X[n+c*R],1,&X[n+c*R],1) * den; }
            }
        }
        free(x1); free(xni);
    }
    else if (iscolmajor && dim==1)
    {
        double *x1, *xni;
        if (!(x1=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)C,&o,0,x1,1); cblas_dcopy((int)C,&ni,0,xni,1);
        for (size_t h=0, n=0, n2=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, n+=RC, n2+=R)
            {
                cblas_dgemv(CblasColMajor,CblasNoTrans,(int)R,(int)C,1.0,&X[n],(int)R,xni,1,0.0,&Y[n2],1);
                cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,1,-1.0,&Y[n2],(int)R,x1,1,1.0,&X[n],(int)R);
                for (size_t r=0; r<R; r++) { Y[n2+r] = cblas_ddot((int)C,&X[n+r],(int)R,&X[n+r],(int)R) * den; }
            }
        }
        free(x1); free(xni);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        double *x1, *xni;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N1,&o,0,x1,1); cblas_dcopy(N1,&ni,0,xni,1);
        for (size_t r=0, n=0, n2=0; r<R; r++, n+=C*SH, n2+=C)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0,&X[n],(int)S,xni,1,0.0,&Y[n2],1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)C,(int)S,1,-1.0,&Y[n2],1,x1,(int)S,1.0,&X[n],(int)S);
            for (size_t c=0; c<C; c++) { Y[n2+c] = cblas_ddot(N1,&X[n+c*S],1,&X[n+c*S],1) * den; }
        }
        free(x1); free(xni);
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (int l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
        {
            for (int m=0; m<M; m++, n+=J, n2++)
            {
                Y[n2] = cblas_ddot(N1,&X[n],K,&ni,0);
                cblas_daxpy(N1,-Y[n2],&o,0,&X[n],K);
                Y[n2] = cblas_ddot(N1,&X[n],K,&X[n],K) * den;
            }
        }
    }
    
    return 0;
}


int var_c (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    if (R<1) { fprintf(stderr,"error in var_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in var_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in var_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in var_c: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = (biased) ? 1.0f/N1 : 1.0f/(N1-1);

    if (N1==1)
    {
        const float z = 0.0f;
        cblas_scopy(N,&z,0,Y,1);
    }
    else if (N1==N)
    {
        const float ni[2] = {-1.0f/N1,0.0f};
        float mn[2] = {0.0f,0.0f};
        for (int n=0; n<2*N1; n+=2) { mn[0] += X[n]; mn[1] += X[n+1]; }
        cblas_caxpy(N1,mn,ni,0,X,1);
        cblas_cdotc_sub(N1,X,1,X,1,(_Complex float *)mn);
        Y[0] = mn[0] * den;
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const float ni[2] = {-1.0f/N1,0.0f};
        float mn[2] = {0.0f,0.0f};
        float *x1;
        if (!(x1=(float *)malloc((size_t)(2*N1)*sizeof(float)))) { fprintf(stderr,"error in var_c: problem with malloc. "); perror("malloc"); return 1; }
        for (int l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
        {
            for (int m=0; m<M; m++, n+=2*J, n2++)
            {
                cblas_ccopy(N1,&X[n],K,x1,1);
                mn[0] = mn[1] = 0.0f;
                for (int n1=0; n1<2*N1; n1+=2) { mn[0] += x1[n1]; mn[1] += x1[n1+1]; }
                cblas_caxpy(N1,mn,ni,0,x1,1);
                cblas_cdotc_sub(N1,x1,1,x1,1,(_Complex float *)mn);
                Y[n2] = mn[0] * den;
            }
        }
        free(x1);
    }

    return 0;
}


int var_z (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    if (R<1) { fprintf(stderr,"error in var_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in var_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in var_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in var_z: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = (biased) ? 1.0/N1 : 1.0/(N1-1);

    if (N1==1)
    {
        const double z = 0.0;
        cblas_dcopy(N,&z,0,Y,1);
    }
    else if (N1==N)
    {
        const double ni[2] = {-1.0/N1,0.0};
        double mn[2] = {0.0,0.0};
        for (int n=0; n<2*N1; n+=2) { mn[0] += X[n]; mn[1] += X[n+1]; }
        cblas_zaxpy(N1,mn,ni,0,X,1);
        cblas_zdotc_sub(N1,X,1,X,1,(_Complex double *)mn);
        Y[0] = mn[0] * den;
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        const double ni[2] = {-1.0/N1,0.0};
        double mn[2] = {0.0,0.0};
        double *x1;
        if (!(x1=(double *)malloc((size_t)(2*N1)*sizeof(double)))) { fprintf(stderr,"error in var_z: problem with malloc. "); perror("malloc"); return 1; }
        for (int l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
        {
            for (int m=0; m<M; m++, n+=2*J, n2++)
            {
                cblas_zcopy(N1,&X[n],K,x1,1);
                mn[0] = mn[1] = 0.0;
                for (int n1=0; n1<2*N1; n1+=2) { mn[0] += x1[n1]; mn[1] += x1[n1+1]; }
                cblas_zaxpy(N1,mn,ni,0,x1,1);
                cblas_zdotc_sub(N1,x1,1,x1,1,(_Complex double *)mn);
                Y[n2] = mn[0] * den;
            }
        }
        free(x1);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
