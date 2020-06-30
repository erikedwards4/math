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

int var_s (float *Y, float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased);
int var_d (double *Y, double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased);
int var_c (float *Y, float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased);
int var_z (double *Y, double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased);


int var_s (float *Y, float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased)
{
    const float z = 0.0f, o = 1.0f;
    float *x1, *xni;
    int r, c, s, h, l, m, n = 0, n2 = 0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (R<1) { fprintf(stderr,"error in var_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in var_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in var_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in var_s: H (num hyperslices X) must be positive\n"); return 1; }

    //Set consts
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ni = 1.0f / N1;
    const float den = (biased) ? 1.0f/N1 : 1.0f/(N1-1);

    //Var
    if (N1==1) { cblas_scopy(N,&z,0,Y,1); }
    else if (N1==N)
    {
        Y[0] = 0.0f;
        while (n<N1) { Y[0] += X[n]; n++; }
        cblas_saxpy(N1,-Y[0],&ni,0,X,1);
        Y[0] = cblas_sdot(N1,X,1,X,1) * den;
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : N1;
        const int ldb = (iscolmajor) ? RC/N1 : 1;
        const int ldc = (iscolmajor) ? R : C;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N1,&ni,0,x1,1);
        cblas_sgemv(Ord,Tr,R,C,1.0f,X,ldc,x1,1,0.0f,Y,1);
        cblas_scopy(N1,&o,0,x1,1);
        if (dim==0)
        {
            cblas_sgemm(Ord,Tr,Tr,R,C,1,-1.0f,x1,lda,Y,ldb,1.0f,X,ldc);
            if (iscolmajor)
            {
                for (c=0; c<C; c++) { Y[c] = cblas_sdot(N1,&X[c*R],1,&X[c*R],1) * den; }
            }
            else
            {
                for (n=0; n<RC; n++) { X[n] *= X[n]; }
                cblas_sgemv(Ord,Tr,R,C,den,X,ldc,x1,1,0.0f,Y,1);
            }
        }
        else
        {
            cblas_sgemm(Ord,Tr,Tr,R,C,1,-1.0f,Y,ldb,x1,lda,1.0f,X,ldc);
            if (iscolmajor)
            {
                for (n=0; n<RC; n++) { X[n] *= X[n]; }
                cblas_sgemv(Ord,Tr,R,C,den,X,ldc,x1,1,0.0f,Y,1);
            }
            else
            {
                for (r=0; r<R; r++) { Y[r] = cblas_sdot(N1,&X[r*C],1,&X[r*C],1) * den; }
            }
        }
        free(x1);
    }
    else if (iscolmajor && dim==0)
    {
        if (!(x1=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(R,&o,0,x1,1); cblas_scopy(R,&ni,0,xni,1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++, n+=RC, n2+=C)
            {
                cblas_sgemv(CblasColMajor,CblasTrans,R,C,1.0f,&X[n],R,xni,1,0.0f,&Y[n2],1);
                cblas_sgemm(CblasColMajor,CblasTrans,CblasTrans,R,C,1,-1.0f,x1,1,&Y[n2],C,1.0f,&X[n],R);
                for (c=0; c<C; c++) { Y[n2+c] = cblas_sdot(R,&X[n+c*R],1,&X[n+c*R],1) * den; }
            }
        }
        free(x1); free(xni);
    }
    else if (iscolmajor && dim==1)
    {
        if (!(x1=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(C,&o,0,x1,1); cblas_scopy(C,&ni,0,xni,1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++, n+=RC, n2+=R)
            {
                cblas_sgemv(CblasColMajor,CblasNoTrans,R,C,1.0f,&X[n],R,xni,1,0.0f,&Y[n2],1);
                cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,R,C,1,-1.0f,&Y[n2],R,x1,1,1.0f,&X[n],R);
                for (r=0; r<R; r++) { Y[n2+r] = cblas_sdot(C,&X[n+r],R,&X[n+r],R) * den; }
            }
        }
        free(x1); free(xni);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N1,&o,0,x1,1); cblas_scopy(N1,&ni,0,xni,1);
        for (r=0; r<R; r++, n+=C*SH, n2+=C)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,C,S,1.0f,&X[n],S,xni,1,0.0f,&Y[n2],1);
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,C,S,1,-1.0f,&Y[n2],1,x1,S,1.0f,&X[n],S);
            for (c=0; c<C; c++) { Y[n2+c] = cblas_sdot(N1,&X[n+c*S],1,&X[n+c*S],1) * den; }
        }
        free(x1); free(xni);
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (l=0; l<L; l++)
        {
            for (m=0, n=l*M*N1; m<M; m++, n+=J, n2++)
            {
                Y[n2] = cblas_sdot(N1,&X[n],K,&ni,0);
                cblas_saxpy(N1,-Y[n2],&o,0,&X[n],K);
                Y[n2] = cblas_sdot(N1,&X[n],K,&X[n],K) * den;
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int var_d (double *Y, double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased)
{
    const double z = 0.0, o = 1.0;
    double *x1, *xni;
    int r, c, s, h, l, m, n = 0, n2 = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in var_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in var_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in var_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in var_d: H (num hyperslices X) must be positive\n"); return 1; }

    //Set consts
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni = 1.0 / N1;
    const double den = (biased) ? 1.0/N1 : 1.0/(N1-1);

    //Var
    if (N1==1) { cblas_dcopy(N,&z,0,Y,1); }
    else if (N1==N)
    {
        Y[0] = 0.0;
        while (n<N1) { Y[0] += X[n]; n++; }
        cblas_daxpy(N1,-Y[0],&ni,0,X,1);
        Y[0] = cblas_ddot(N1,X,1,X,1) * den;
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : N1;
        const int ldb = (iscolmajor) ? RC/N1 : 1;
        const int ldc = (iscolmajor) ? R : C;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N1,&ni,0,x1,1);
        cblas_dgemv(Ord,Tr,R,C,1.0,X,ldc,x1,1,0.0,Y,1);
        cblas_dcopy(N1,&o,0,x1,1);
        if (dim==0)
        {
            cblas_dgemm(Ord,Tr,Tr,R,C,1,-1.0,x1,lda,Y,ldb,1.0,X,ldc);
            if (iscolmajor)
            {
                for (c=0; c<C; c++) { Y[c] = cblas_ddot(N1,&X[c*R],1,&X[c*R],1) * den; }
            }
            else
            {
                for (n=0; n<RC; n++) { X[n] *= X[n]; }
                cblas_dgemv(Ord,Tr,R,C,den,X,ldc,x1,1,0.0,Y,1);
            }
        }
        else
        {
            cblas_dgemm(Ord,Tr,Tr,R,C,1,-1.0,Y,ldb,x1,lda,1.0,X,ldc);
            if (iscolmajor)
            {
                for (n=0; n<RC; n++) { X[n] *= X[n]; }
                cblas_dgemv(Ord,Tr,R,C,den,X,ldc,x1,1,0.0,Y,1);
            }
            else
            {
                for (r=0; r<R; r++) { Y[r] = cblas_ddot(N1,&X[r*C],1,&X[r*C],1) * den; }
            }
        }
        free(x1);
    }
    else if (iscolmajor && dim==0)
    {
        if (!(x1=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(R,&o,0,x1,1); cblas_dcopy(R,&ni,0,xni,1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++, n+=RC, n2+=C)
            {
                cblas_dgemv(CblasColMajor,CblasTrans,R,C,1.0,&X[n],R,xni,1,0.0,&Y[n2],1);
                cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,R,C,1,-1.0,x1,1,&Y[n2],C,1.0,&X[n],R);
                for (c=0; c<C; c++) { Y[n2+c] = cblas_ddot(R,&X[n+c*R],1,&X[n+c*R],1) * den; }
            }
        }
        free(x1); free(xni);
    }
    else if (iscolmajor && dim==1)
    {
        if (!(x1=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(C,&o,0,x1,1); cblas_dcopy(C,&ni,0,xni,1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++, n+=RC, n2+=R)
            {
                cblas_dgemv(CblasColMajor,CblasNoTrans,R,C,1.0,&X[n],R,xni,1,0.0,&Y[n2],1);
                cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,R,C,1,-1.0,&Y[n2],R,x1,1,1.0,&X[n],R);
                for (r=0; r<R; r++) { Y[n2+r] = cblas_ddot(C,&X[n+r],R,&X[n+r],R) * den; }
            }
        }
        free(x1); free(xni);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N1,&o,0,x1,1); cblas_dcopy(N1,&ni,0,xni,1);
        for (r=0; r<R; r++, n+=C*SH, n2+=C)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,C,S,1.0,&X[n],S,xni,1,0.0,&Y[n2],1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,C,S,1,-1.0,&Y[n2],1,x1,S,1.0,&X[n],S);
            for (c=0; c<C; c++) { Y[n2+c] = cblas_ddot(N1,&X[n+c*S],1,&X[n+c*S],1) * den; }
        }
        free(x1); free(xni);
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (l=0; l<L; l++)
        {
            for (m=0, n=l*M*N1; m<M; m++, n+=J, n2++)
            {
                Y[n2] = cblas_ddot(N1,&X[n],K,&ni,0);
                cblas_daxpy(N1,-Y[n2],&o,0,&X[n],K);
                Y[n2] = cblas_ddot(N1,&X[n],K,&X[n],K) * den;
            }
        }
    }

    return 0;
}


int var_c (float *Y, float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased)
{
    const float z[2] =  {0.0f,0.0f};
    float mn[2] = {0.0f,0.0f};
    float *x1;
    int l, m, n = 0, n1, n2 = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in var_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in var_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in var_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in var_c: H (num hyperslices X) must be positive\n"); return 1; }

    //Set consts
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ni[2] = {-1.0f/N1,0.0f};
    const float den = (biased) ? 1.0f/N1 : 1.0f/(N1-1);

    //Var
    if (N1==1) { cblas_scopy(N,&z[0],0,Y,1); }
    else if (N1==N)
    {
        while (n<2*N1) { mn[0] += X[n]; mn[1] += X[n+1]; n+=2; }
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
        if (!(x1=(float *)malloc((size_t)(2*N1)*sizeof(float)))) { fprintf(stderr,"error in var_c: problem with malloc. "); perror("malloc"); return 1; }
        for (l=0; l<L; l++)
        {
            for (m=0, n=l*M*N1; m<M; m++, n+=J, n2++)
            {
                cblas_ccopy(N1,&X[2*n],K,x1,1);
                mn[0] = mn[1] = 0.0f;
                for (n1=0; n1<2*N1; n1+=2) { mn[0] += x1[n1]; mn[1] += x1[n1+1]; }
                cblas_caxpy(N1,mn,ni,0,x1,1);
                cblas_cdotc_sub(N1,x1,1,x1,1,(_Complex float *)mn);
                Y[n2] = mn[0] * den;
            }
        }
    }

    return 0;
}


int var_z (double *Y, double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased)
{
    const double z[2] =  {0.0,0.0};
    double mn[2] = {0.0,0.0};
    double *x1;
    int l, m, n = 0, n1, n2 = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in var_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in var_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in var_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in var_z: H (num hyperslices X) must be positive\n"); return 1; }

    //Set consts
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni[2] = {-1.0/N1,0.0};
    const double den = (biased) ? 1.0/N1 : 1.0/(N1-1);

    //Var
    if (N1==1) { cblas_dcopy(N,&z[0],0,Y,1); }
    else if (N1==N)
    {
        while (n<2*N1) { mn[0] += X[n]; mn[1] += X[n+1]; n+=2; }
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
        if (!(x1=(double *)malloc((size_t)(2*N1)*sizeof(double)))) { fprintf(stderr,"error in var_z: problem with malloc. "); perror("malloc"); return 1; }
        for (l=0; l<L; l++)
        {
            for (m=0, n=l*M*N1; m<M; m++, n+=J, n2++)
            {
                cblas_zcopy(N1,&X[2*n],K,x1,1);
                mn[0] = mn[1] = 0.0;
                for (n1=0; n1<2*N1; n1+=2) { mn[0] += x1[n1]; mn[1] += x1[n1+1]; }
                cblas_zaxpy(N1,mn,ni,0,x1,1);
                cblas_zdotc_sub(N1,x1,1,x1,1,(_Complex double *)mn);
                Y[n2] = mn[0] * den;
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
