//Gets sum of each row or col of X according to dim.
//For complex case, real and imag parts calculated separately.

//For 2D case, this is definitely faster using cblas_?gemv.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sum_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor);
int sum_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor);
int sum_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor);
int sum_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor);


int sum_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor)
{
    const float o = 1.0f;
    float *x1;
    int r, s, h, l, m, n = 0, n2 = 0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (R<1) { fprintf(stderr,"error in sum_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in sum_s: C (ncols X) must be positive\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sum_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sum_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_scopy(N,X,1,Y,1); }
    else if (N1==N) { Y[0] = 0.0f; while (n<N) { Y[0] += X[n]; n++; } }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? R : C;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in sum_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N1,&o,0,x1,1);
        cblas_sgemv(Ord,Tr,R,C,1.0f,X,lda,x1,1,0.0f,Y,1);
        //for (c=0; c<C; c++) { Y[c] = cblas_sdot(R,&X[c*R],1,&o,0); }
        //for (n=0, c=0; n<RC; n++, c++) { if (c==C) { c=0; } Y[c] += X[n]; }
        free(x1);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in sum_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N1,&o,0,x1,1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++, n+=RC, n2+=RC/N1)
            {
                cblas_sgemv(CblasColMajor,Tr,R,C,1.0f,&X[n],R,x1,1,0.0f,&Y[n2],1);
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in sum_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N1,&o,0,x1,1);
        for (r=0; r<R; r++, n+=C*SH, n2+=C)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,C,S,1.0f,&X[n],S,x1,1,0.0f,&Y[n2],1);
        }
        free(x1);
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
                Y[n2] = cblas_sdot(N1,&X[n],K,&o,0);
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int sum_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor)
{
    const double o = 1.0;
    double *x1;
    int r, s, h, l, m, n = 0, n2 = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in sum_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in sum_d: C (ncols X) must be positive\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sum_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sum_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_dcopy(N,X,1,Y,1); }
    else if (N1==N) { Y[0] = 0.0; while (n<N) { Y[0] += X[n]; n++; } }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? R : C;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in sum_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N1,&o,0,x1,1);
        cblas_dgemv(Ord,Tr,R,C,1.0,X,lda,x1,1,0.0,Y,1);
        free(x1);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in sum_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N1,&o,0,x1,1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++, n+=RC, n2+=RC/N1)
            {
                cblas_dgemv(CblasColMajor,Tr,R,C,1.0,&X[n],R,x1,1,0.0,&Y[n2],1);
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in sum_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N1,&o,0,x1,1);
        for (r=0; r<R; r++, n+=C*SH, n2+=C)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,C,S,1.0,&X[n],S,x1,1,0.0,&Y[n2],1);
        }
        free(x1);
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
                Y[n2] = cblas_ddot(N1,&X[n],K,&o,0);
            }
        }
    }
    
    return 0;
}


int sum_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor)
{
    const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
    float *x1;
    int r, s, h, l, m, n = 0, n2 = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in sum_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in sum_c: C (ncols X) must be positive\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sum_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sum_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_ccopy(N,X,1,Y,1); }
    else if (N1==N)
    {
        Y[0] = Y[1] = 0.0f;
        while (n<2*N1) { Y[0] += X[n]; Y[1] += X[n+1]; n+=2; }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? R : C;
        if (!(x1=(float *)malloc((size_t)(2*N1)*sizeof(float)))) { fprintf(stderr,"error in sum_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy(N1,o,0,x1,1);
        cblas_cgemv(Ord,Tr,R,C,o,X,lda,x1,1,z,Y,1);
        free(x1);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        if (!(x1=(float *)malloc((size_t)(2*N1)*sizeof(float)))) { fprintf(stderr,"error in sum_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy(N1,o,0,x1,1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++, n+=2*RC, n2+=2*RC/N1)
            {
                cblas_cgemv(CblasColMajor,Tr,R,C,o,&X[n],R,x1,1,z,&Y[n2],1);
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        if (!(x1=(float *)malloc((size_t)(2*N1)*sizeof(float)))) { fprintf(stderr,"error in sum_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy(N1,o,0,x1,1);
        for (r=0; r<R; r++, n+=2*C*SH, n2+=2*C)
        {
            cblas_cgemv(CblasRowMajor,CblasNoTrans,C,S,o,&X[n],S,x1,1,z,&Y[n2],1);
        }
        free(x1);
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (l=0; l<L; l++)
        {
            for (m=0, n=2*l*M*N1; m<M; m++, n+=2*J, n2+=2)
            {
                cblas_cdotu_sub(N1,&X[n],K,o,0,(openblas_complex_float *)(&Y[n2]));
            }
        }
    }
    
    return 0;
}


int sum_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor)
{
    const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
    double *x1;
    int r, s, h, l, m, n = 0, n2 = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in sum_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in sum_z: C (ncols X) must be positive\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sum_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sum_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_zcopy(N,X,1,Y,1); }
    else if (N1==N)
    {
        Y[0] = Y[1] = 0.0;
        while (n<2*N1) { Y[0] += X[n]; Y[1] += X[n+1]; n+=2; }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? R : C;
        if (!(x1=(double *)malloc((size_t)(2*N1)*sizeof(double)))) { fprintf(stderr,"error in sum_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy(N1,o,0,x1,1);
        cblas_zgemv(Ord,Tr,R,C,o,X,lda,x1,1,z,Y,1);
        free(x1);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        if (!(x1=(double *)malloc((size_t)(2*N1)*sizeof(double)))) { fprintf(stderr,"error in sum_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy(N1,o,0,x1,1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++, n+=2*RC, n2+=2*RC/N1)
            {
                cblas_zgemv(CblasColMajor,Tr,R,C,o,&X[n],R,x1,1,z,&Y[n2],1);
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        if (!(x1=(double *)malloc((size_t)(2*N1)*sizeof(double)))) { fprintf(stderr,"error in sum_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy(N1,o,0,x1,1);
        for (r=0; r<R; r++, n+=2*C*SH, n2+=2*C)
        {
            cblas_zgemv(CblasRowMajor,CblasNoTrans,C,S,o,&X[n],S,x1,1,z,&Y[n2],1);
        }
        free(x1);
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (l=0; l<L; l++)
        {
            for (m=0, n=2*l*M*N1; m<M; m++, n+=2*J, n2+=2)
            {
                cblas_zdotu_sub(N1,&X[n],K,o,0,(openblas_complex_double *)(&Y[n2]));
            }
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
