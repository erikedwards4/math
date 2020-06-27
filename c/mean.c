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

int mean_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor);
int mean_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor);
int mean_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor);
int mean_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor);


int mean_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor)
{
    float *x1;
    int r, s, h, l, m, n1 = 0, n2 = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in mean_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mean_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in mean_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in mean_s: H (num hyperslices X) must be positive\n"); return 1; }

    //Set ints
    const int RC = R*C, SH = S*H;
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = RC*SH/(M*N);
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const float ni = 1.0f / N;

    if (N==1) { cblas_scopy(RC*SH,X,1,Y,1); }
    else if (N==RC*SH)
    {
        Y[0] = 0.0f;
        while (n1<N) { Y[0] += X[n1]; n1++; }
        Y[0] *= ni;
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? R : C;
        if (!(x1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in mean_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N,&ni,0,x1,1);
        cblas_sgemv(Ord,Tr,R,C,1.0f,X,lda,x1,1,0.0f,Y,1);
        free(x1);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        if (!(x1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in mean_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N,&ni,0,x1,1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++, n1+=RC, n2+=RC/N)
            {
                cblas_sgemv(CblasColMajor,Tr,R,C,1.0f,&X[n1],R,x1,1,0.0f,&Y[n2],1);
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        if (!(x1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in mean_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N,&ni,0,x1,1);
        for (r=0; r<R; r++, n1+=C*SH, n2+=C)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,C,S,1.0f,&X[n1],S,x1,1,0.0f,&Y[n2],1);
        }
        free(x1);
    }
    else
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n1=l*M*N; m<M; m++, n1+=J, n2++)
            {
                Y[n2] = cblas_sdot(N,&X[n1],I,&ni,0);
            }
        }
    }
    
    return 0;
}


int mean_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor)
{
    double *x1;
    int r, s, h, l, m, n1 = 0, n2 = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in mean_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mean_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in mean_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in mean_d: H (num hyperslices X) must be nonnepositivegative\n"); return 1; }

    //Set ints
    const int RC = R*C, SH = S*H;
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = RC*SH/(M*N);
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const double ni = 1.0 / N;

    if (N==1) { cblas_dcopy(RC*SH,X,1,Y,1); }
    else if (N==RC*SH)
    {
        Y[0] = 0.0;
        while (n1<N) { Y[0] += X[n1]; n1++; }
        Y[0] *= ni;
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? R : C;
        if (!(x1=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in mean_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N,&ni,0,x1,1);
        cblas_dgemv(Ord,Tr,R,C,1.0,X,lda,x1,1,0.0,Y,1);
        free(x1);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        if (!(x1=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in mean_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N,&ni,0,x1,1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++, n1+=RC, n2+=RC/N)
            {
                cblas_dgemv(CblasColMajor,Tr,R,C,1.0,&X[n1],R,x1,1,0.0,&Y[n2],1);
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        if (!(x1=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in mean_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N,&ni,0,x1,1);
        for (r=0; r<R; r++, n1+=C*SH, n2+=C)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,C,S,1.0,&X[n1],S,x1,1,0.0,&Y[n2],1);
        }
        free(x1);
    }
    else
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n1=l*M*N; m<M; m++, n1+=J, n2++)
            {
                Y[n2] = cblas_ddot(N,&X[n1],I,&ni,0);
            }
        }
    }
    
    return 0;
}


int mean_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor)
{
    const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
    float *x1;
    int r, s, h, l, m, n1 = 0, n2 = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in mean_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mean_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in mean_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in mean_c: H (num hyperslices X) must be positive\n"); return 1; }

    //Set ints
    const int RC = R*C, SH = S*H;
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = RC*SH/(M*N);
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const float ni[2] = {1.0f/N,0.0f};

    if (N==1) {}
    else if (N==RC*SH)
    {
        Y[0] = Y[1] = 0.0f;
        while (n1<2*N) { Y[0] += X[n1]; Y[1] += X[n1+1]; n1+=2; }
        Y[0] *= ni[0]; Y[1] *= ni[0];
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? R : C;
        if (!(x1=(float *)malloc((size_t)(2*N)*sizeof(float)))) { fprintf(stderr,"error in mean_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy(N,ni,0,x1,1);
        cblas_cgemv(Ord,Tr,R,C,o,X,lda,x1,1,z,Y,1);
        free(x1);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        if (!(x1=(float *)malloc((size_t)(2*N)*sizeof(float)))) { fprintf(stderr,"error in mean_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy(N,ni,0,x1,1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++, n1+=2*RC, n2+=2*RC/N)
            {
                cblas_cgemv(CblasColMajor,Tr,R,C,o,&X[n1],R,x1,1,z,&Y[n2],1);
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        if (!(x1=(float *)malloc((size_t)(2*N)*sizeof(float)))) { fprintf(stderr,"error in mean_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy(N,ni,0,x1,1);
        for (r=0; r<R; r++, n1+=2*C*SH, n2+=2*C)
        {
            cblas_cgemv(CblasRowMajor,CblasNoTrans,C,S,o,&X[n1],S,x1,1,z,&Y[n2],1);
        }
        free(x1);
    }
    else
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n1=2*l*M*N; m<M; m++, n1+=2*J, n2+=2)
            {
                cblas_cdotu_sub(N,&X[n1],I,ni,0,(openblas_complex_float *)(&Y[n2]));
            }
        }
    }
    
    return 0;
}


int mean_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const int dim, const char iscolmajor)
{
    const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
    double *x1;
    int r, s, h, l, m, n1 = 0, n2 = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in mean_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mean_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in mean_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in mean_z: H (num hyperslices X) must be positive\n"); return 1; }

    //Set ints
    const int RC = R*C, SH = S*H;
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = RC*SH/(M*N);
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const double ni[2] = {1.0/N,0.0};

    if (N==1) {}
    else if (N==RC*SH)
    {
        Y[0] = Y[1] = 0.0;
        while (n1<2*N) { Y[0] += X[n1]; Y[1] += X[n1+1]; n1+=2; }
        Y[0] *= ni[0]; Y[1] *= ni[0];
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? R : C;
        if (!(x1=(double *)malloc((size_t)(2*N)*sizeof(double)))) { fprintf(stderr,"error in mean_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy(N,ni,0,x1,1);
        cblas_zgemv(Ord,Tr,R,C,o,X,lda,x1,1,z,Y,1);
        free(x1);
    }
    else if (iscolmajor && dim<2)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        if (!(x1=(double *)malloc((size_t)(2*N)*sizeof(double)))) { fprintf(stderr,"error in mean_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy(N,ni,0,x1,1);
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++, n1+=2*RC, n2+=2*RC/N)
            {
                cblas_zgemv(CblasColMajor,Tr,R,C,o,&X[n1],R,x1,1,z,&Y[n2],1);
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        if (!(x1=(double *)malloc((size_t)(2*N)*sizeof(double)))) { fprintf(stderr,"error in mean_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy(N,ni,0,x1,1);
        for (r=0; r<R; r++, n1+=2*C*SH, n2+=2*C)
        {
            cblas_zgemv(CblasRowMajor,CblasNoTrans,C,S,o,&X[n1],S,x1,1,z,&Y[n2],1);
        }
        free(x1);
    }
    else
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n1=2*l*M*N; m<M; m++, n1+=2*J, n2+=2)
            {
                cblas_zdotu_sub(N,&X[n1],I,ni,0,(openblas_complex_double *)(&Y[n2]));
            }
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
