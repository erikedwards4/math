//Gets skewness of each row or col of X according to dim.
//For complex case, output is complex.
//I follow the Octave convention for complex skewness (but see literature for other ideas later).
//Complex skewness only supported for vectors currently.
//This works in place.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int skewness_s (float *Y, float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased);
int skewness_d (double *Y, double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased);
int skewness_c (float *Y, float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased);
int skewness_z (double *Y, double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased);


int skewness_s (float *Y, float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased)
{
    const float o = 1.0f;
    int l, m, r, c, n1 = 0, n2 = 0;
    float x, *x1, *xni;
    struct timespec tic, toc;
    clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (R<1) { fprintf(stderr,"error in skewness_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in skewness_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in skewness_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in skewness_s: H (num hyperslices X) must be positive\n"); return 1; }

    //Set consts
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int N2 = N/N1;
    const float ni = 1.0f / N;
    const float w = (biased) ? sqrtf(N1) : N*sqrtf(N1-1)/(N1-2);

    //Skewness
    if (N1<2) { fprintf(stderr,"error in skewness_s: N must be > 1\n"); return 1; }
    else if (N1==N)
    {
        float sm2 = 0.0f, sm3 = 0.0f;
        Y[0] = 0.0f;
        while (n1<N1) { Y[0] += X[n1]; n1++; }
        cblas_saxpy(N,-Y[0],&ni,0,X,1);
        while (n2<N1) { x = X[n2]*X[n2]; sm2 += x; sm3 += x*X[n2]; n2++; }
        Y[0] = w * sm3 / (sm2*sqrtf(sm2));
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : N;
        const int ldb = (iscolmajor) ? RC/N : 1;
        const int ldc = (iscolmajor) ? R : C;
        float *sm2, *sm3;
        if (!(x1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N,&ni,0,x1,1);
        cblas_sgemv(Ord,Tr,R,C,1.0f,X,ldc,x1,1,0.0f,Y,1);
        cblas_scopy(N,&o,0,x1,1);
        if (dim==0) { cblas_sgemm(Ord,Tr,Tr,R,C,1,-1.0f,x1,lda,Y,ldb,1.0f,X,ldc); }
        else { cblas_sgemm(Ord,Tr,Tr,R,C,1,-1.0f,Y,ldb,x1,lda,1.0f,X,ldc); }
        if (dim==0)
        {
            if (iscolmajor)
            {
                for (c=0; c<C; c++) { Y[c] = cblas_sdot(N,&X[c*R],1,&X[c*R],1); }
                while (n1<RC) { X[n1] *= X[n1]*X[n1]; n1++; }
                for (c=0; c<C; c++) { Y[c] = w * cblas_sdot(N,&X[c*R],1,x1,1) / (Y[c]*sqrtf(Y[c])); }
                //while (n1<RC) { c = n1%R; x = X[n1]*X[n1]; sm2[c] += x; sm3[c] += x*X[n1]; n1++; }
                //for (c=0; c<C; c++) { Y[c] = w * sm3[c] / (sm2[c]*sqrtf(sm2[c])); }
            }
            else
            {
                //for (c=0; c<C; c++) { Y[c] = cblas_sdot(N,&X[c],C,&X[c],C); }
                //while (n1<RC) { X[n1] *= X[n1]*X[n1]; n1++; }
                //while (n1<RC) { x = X[n1]*X[n1]; Y[n1%C] += x; X[n1] *= x; n1++; }
                //for (c=0; c<C; c++) { Y[c] = w * cblas_sdot(N,&X[c],C,x1,1) / (Y[c]*sqrtf(Y[c])); }
                if (!(sm2=(float *)calloc((size_t)N,sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with calloc. "); perror("calloc"); return 1; }
                if (!(sm3=(float *)calloc((size_t)N,sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with calloc. "); perror("calloc"); return 1; }
                while (n1<RC) { c = n1%R; x = X[n1]*X[n1]; sm2[c] += x; sm3[c] += x*X[n1]; n1++; }
                for (c=0; c<C; c++) { Y[c] = w * sm3[c] / (sm2[c]*sqrtf(sm2[c])); }
                free(sm2); free(sm3);
            }
        }
        else
        {
            if (iscolmajor)
            {
                if (!(sm2=(float *)calloc((size_t)N,sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with calloc. "); perror("calloc"); return 1; }
                if (!(sm3=(float *)calloc((size_t)N,sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with calloc. "); perror("calloc"); return 1; }
                while (n1<RC) { r = n1%C; x = X[n1]*X[n1]; sm2[r] += x; sm3[r] += x*X[n1]; n1++; }
                for (r=0; r<R; r++) { Y[r] = w * sm3[r] / (sm2[r]*sqrtf(sm2[r])); }
                free(sm2); free(sm3);
            }
            else
            {
                for (r=0; r<R; r++) { Y[r] = cblas_sdot(N,&X[r*C],1,&X[r*C],1); }
                while (n1<RC) { X[n1] *= X[n1]*X[n1]; n1++; }
                for (r=0; r<R; r++) { Y[r] = w * cblas_sdot(N,&X[r*C],1,x1,1) / (Y[r]*sqrtf(Y[r])); }
            }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        float *sm2, *sm3;
        if (!(sm2=(float *)calloc((size_t)N2,sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with calloc. "); perror("calloc"); return 1; }
        if (!(sm3=(float *)calloc((size_t)N2,sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with calloc. "); perror("calloc"); return 1; }
        if (!(x1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N,&o,0,x1,1); cblas_scopy(N,&ni,0,xni,1);
        for (r=0; r<R; r++, n1+=C*SH, n2+=C)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,C,S,1.0f,&X[n1],S,xni,1,0.0f,&Y[n2],1);
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,C,S,1,-1.0f,&Y[n2],1,x1,S,1.0f,&X[n1],S);
        }
        for (n1=0; n1<N; n1++) { n2 = n1/N; x = X[n1]*X[n1]; sm2[n2] += x; sm3[n2] += x*X[n1]; }
        for (n2=0; n2<RC; n2++) { Y[n2] = w * sm3[n2] / (sm2[n2]*sqrtf(sm2[n2])); }
        free(x1); free(xni); free(sm2); free(sm3);
    }
    else
    {
        fprintf(stderr,"else\n");
        float *sm2, *sm3;
        if (!(sm2=(float *)calloc((size_t)N2,sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with calloc. "); perror("calloc"); return 1; }
        if (!(sm3=(float *)calloc((size_t)N2,sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with calloc. "); perror("calloc"); return 1; }
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        for (l=0; l<L; l++)
        {
            for (m=0, n1=l*M*N; m<M; m++, n1+=J, n2++)
            {
                Y[n2] = cblas_sdot(N,&X[n1],K,&ni,0);
                cblas_saxpy(N,-Y[n2],&o,0,&X[n1],K);
                Y[n2] = cblas_sdot(N,&X[n1],K,&X[n1],K);
            }
        }
        //for (n1=0; n1<N; n1++) { n2 = n1/N; x = X[n1]*X[n1]; sm2[n2] += x; sm3[n2] += x*X[n1]; }
        for (n1=0; n1<N; n1++) { X[n1] *= X[n1]*X[n1]; }
        for (l=0; l<L; l++)
        {
            for (m=0, n1=l*M*N; m<M; m++, n1+=J, n2++)
            {
                Y[n2] = w * cblas_sdot(N,&X[n1],K,&X[n1],K) / (Y[n2]*sqrtf(Y[n2]));
            }
        }
    }
    clock_gettime(CLOCK_REALTIME,&toc);
    fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int skewness_d (double *Y, double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased)
{
    const double o = 1.0;
    int r, c, l, m, n, n1 = 0, n2 = 0;
    double *x1, *xni;

    //Checks
    if (R<1) { fprintf(stderr,"error in skewness_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in skewness_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in skewness_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in skewness_d: H (num hyperslices X) must be positive\n"); return 1; }

    //Set consts
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni = 1.0 / N;

    //Skewness
    if (N1<2) { fprintf(stderr,"error in skewness_d: N must be > 1\n"); return 1; }
    else if (N1==N)
    {
        double sm2 = 0.0, sm3 = 0.0;
        Y[0] = 0.0;
        while (n1<N1) { Y[0] += X[n1]; n1++; }
        cblas_daxpy(N,-Y[0],&ni,0,X,1);
        while (n2<N1) { X[n2] *= X[n2]; sm2 += X[n2]; sm3 += X[n2]*X[n2]; n2++; }
        Y[0] = N*sm3/(sm2*sm2);
        if (!biased) { Y[0] =  3.0 + (Y[0]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : N;
        const int ldb = (iscolmajor) ? RC/N : 1;
        const int ldc = (iscolmajor) ? R : C;
        if (!(x1=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N,&ni,0,x1,1);
        cblas_dgemv(Ord,Tr,R,C,1.0,X,ldc,x1,1,0.0,Y,1);
        cblas_dcopy(N,&o,0,x1,1);
        if (dim==0) { cblas_dgemm(Ord,Tr,Tr,R,C,1,-1.0,x1,lda,Y,ldb,1.0,X,ldc); }
        else { cblas_dgemm(Ord,Tr,Tr,R,C,1,-1.0,Y,ldb,x1,lda,1.0,X,ldc); }
        if (dim==0)
        {
            if (iscolmajor)
            {
                while (n1<RC) { X[n1] *= X[n1]; n1++; }
                for (c=0; c<C; c++)
                {
                    Y[c] = cblas_ddot(N,&X[c*R],1,x1,1);
                    Y[c] = N * cblas_ddot(N,&X[c*R],1,&X[c*R],1) / (Y[c]*Y[c]);
                }
            }
            else
            {
                double *sm2, *sm3;
                if (!(sm2=(double *)calloc((size_t)C,sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with calloc. "); perror("calloc"); return 1; }
                if (!(sm3=(double *)calloc((size_t)C,sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with calloc. "); perror("calloc"); return 1; }
                while (n1<RC) { c = n1%C; X[n1] *= X[n1]; sm2[c] += X[n1]; sm3[c] += X[n1]*X[n1]; n1++; }
                for (c=0; c<C; c++) { Y[c] = N*sm3[c]/(sm2[c]*sm2[c]); }
                free(sm2); free(sm3);
            }
        }
        else
        {
            if (iscolmajor)
            {
                double *sm2, *sm3;
                if (!(sm2=(double *)calloc((size_t)R,sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with calloc. "); perror("calloc"); return 1; }
                if (!(sm3=(double *)calloc((size_t)R,sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with calloc. "); perror("calloc"); return 1; }
                while (n1<RC) { r = n1%R; X[n1] *= X[n1]; sm2[r] += X[n1]; sm3[r] += X[n1]*X[n1]; n1++; }
                for (r=0; r<R; r++) { Y[r] = N*sm3[r]/(sm2[r]*sm2[r]); }
                free(sm2); free(sm3);
            }
            else
            {
                while (n1<RC) { X[n1] *= X[n1]; n1++; }
                for (r=0; r<R; r++)
                {
                    Y[r] = cblas_ddot(N,&X[r*C],1,x1,1);
                    Y[r] = N * cblas_ddot(N,&X[r*C],1,&X[r*C],1) / (Y[r]*Y[r]);
                }
            }
        }
        if (!biased)
        {
            for (n=0; n<RC/N; n++) { Y[n] =  3.0 + (Y[n]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        if (!(x1=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N,&o,0,x1,1); cblas_dcopy(N,&ni,0,xni,1);
        for (r=0; r<R; r++, n1+=C*SH, n2+=C)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,C,S,1.0,&X[n1],S,xni,1,0.0,&Y[n2],1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,C,S,1,-1.0,&Y[n2],1,x1,S,1.0,&X[n1],S);
        }
        for (n=0; n<N; n++) { X[n] *= X[n]; }
        for (r=0, n1=0, n2=0; r<R; r++, n1+=C*SH, n2+=C)
        {
            for (c=0; c<C; c++)
            {
                Y[n2+c] = cblas_ddot(N,&X[n1+c*S],1,x1,1);
                Y[n2+c] = N * cblas_ddot(N,&X[n1+c*S],1,&X[n1+c*S],1) / (Y[n2+c]*Y[n2+c]);
            }
        }
        if (!biased)
        {
            for (n=0; n<RC; n++) { Y[n] =  3.0 + (Y[n]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1); free(xni);
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        if (!(x1=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N,&o,0,x1,1);
        for (l=0; l<L; l++)
        {
            for (m=0, n1=l*M*N; m<M; m++, n1+=J, n2++)
            {
                cblas_ddot(N,&X[n1],K,&ni,0);
                cblas_daxpy(N,-Y[n2],&o,0,&X[n1],K);
                for (n=0; n<N; n++) { X[n1+n*K] *= X[n1+n*K]; }
                Y[n2] = cblas_ddot(N,&X[n1],K,x1,1);
                Y[n2] = N * cblas_ddot(N,&X[n1],K,&X[n1],K) / (Y[n2]*Y[n2]);
            }
        }
        if (!biased)
        {
            for (n=0; n<N/N; n++) { Y[n] =  3.0 + (Y[n]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1);
    }

    return 0;
}


int skewness_c (float *Y, float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased)
{
    const float z[2] =  {0.0f,0.0f};
    float mn[2] = {0.0f,0.0f}, sm2[2] = {0.0f,0.0f}, sm3[2] = {0.0f,0.0f};
    float tmp, *x1;
    int n = 0, n1 = 0;
    //int l, m, n = 0, n1 = 0, n2 = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in skewness_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in skewness_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in skewness_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in skewness_c: H (num hyperslices X) must be positive\n"); return 1; }

    //Set consts
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ni[2] = {-1.0f/N,0.0f};

    //Allocate
    if (!(x1=(float *)malloc((size_t)(2*N1)*sizeof(float)))) { fprintf(stderr,"error in skewness_c: problem with malloc. "); perror("malloc"); return 1; }

    //Skewness
    if (N1==1) { cblas_ccopy(N,z,0,Y,1); }
    else if (N1==N)
    {
        while (n1<2*N1) { mn[0] += X[n1]; mn[1] += X[n1+1]; n1+=2; }
        cblas_caxpy(N,mn,ni,0,X,1);
        cblas_cdotc_sub(N,X,1,X,1,(_Complex float *)sm2);
        for (n=0; n<2*N; n+=2)
        {
            x1[n] = X[n]*X[n] - X[n+1]*X[n+1];
            x1[n+1] = X[n]*X[n+1] + X[n+1]*X[n];
        }
        for (n=0; n<2*N; n+=2)
        {
            tmp = x1[n]*X[n] - x1[n+1]*X[n+1];
            x1[n+1] = x1[n]*X[n+1] + x1[n+1]*X[n];
            x1[n] = tmp;
        }
        cblas_cdotu_sub(N,x1,1,X,1,(_Complex float *)sm3);
        Y[0] = N*sm3[0]/(sm2[0]*sm2[0]); Y[1] = N*sm3[1]/(sm2[0]*sm2[0]);
        if (!biased)
        {
            Y[0] =  3.0f + (Y[0]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3));
            Y[1] *= (N1+1)*(N1-1)/((N1-2)*(N1-3));
        }
    }
    else
    {
        if (iscolmajor) {}
        fprintf(stderr,"error in skewness_c: complex skewness only supported for vectors currently \n"); return 1;
        // for (l=0; l<L; l++)
        // {
        //     for (m=0, n1=l*M*N; m<M; m++, n1+=J, n2+=2)
        //     {
        //         cblas_ccopy(N,&X[2*n1],K,x1,1);
        //         mn[0] = mn[1] = 0.0f;
        //         for (n=0; n<2*N; n+=2) { mn[0] += x1[n]; mn[1] += x1[n+1]; }
        //         cblas_caxpy(N,mn,ni,0,&X[2*n1],K);
        //         cblas_caxpy(N,mn,ni,0,x1,1);
        //         cblas_cdotc_sub(N,x1,1,x1,1,(_Complex float *)sm2);
        //         for (n=0; n<2*N; n+=2)
        //         {
        //             tmp = x1[n]*x1[n] - x1[n+1]*x1[n+1];
        //             x1[n+1] = x1[n]*x1[n+1] + x1[n+1]*x1[n];
        //             x1[n] = tmp;
        //         }
        //         for (n=0; n<2*N; n+=2)
        //         {
        //             tmp = x1[n]*X[2*n1+n*K] - x1[n+1]*X[2*n1+n*K+1];
        //             x1[n+1] = x1[n]*X[2*n1+n*K+1] + x1[n+1]*X[2*n1+n*K];
        //             x1[n] = tmp;
        //         }
        //         cblas_cdotu_sub(N,x1,1,&X[2*n1],K,(_Complex float *)sm3);
        //         Y[n2] = N*sm3[0]/(sm2[0]*sm2[0]); Y[n2+1] = N*sm3[1]/(sm2[0]*sm2[0]);
        //     }
        // }
        // if (!biased)
        // {
        //     for (n=0; n<N/N; n+=2)
        //     {
        //         Y[n] =  3.0f + (Y[n]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3));
        //         Y[n+1] *= (N1+1)*(N1-1)/((N1-2)*(N1-3));
        //     }
        // }
    }

    return 0;
}


int skewness_z (double *Y, double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor, const char biased)
{
    const double z[2] =  {0.0,0.0};
    double mn[2] = {0.0,0.0}, sm2[2] = {0.0,0.0}, sm3[2] = {0.0,0.0};
    double tmp, *x1;
    int n = 0, n1 = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in skewness_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in skewness_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in skewness_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in skewness_z: H (num hyperslices X) must be positive\n"); return 1; }

    //Set consts
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni[2] = {-1.0/N,0.0};

    //Allocate
    if (!(x1=(double *)malloc((size_t)(2*N1)*sizeof(double)))) { fprintf(stderr,"error in skewness_z: problem with malloc. "); perror("malloc"); return 1; }

    //Skewness
    if (N1==1) { cblas_zcopy(N,z,0,Y,1); }
    else if (N1==N)
    {
        while (n1<2*N1) { mn[0] += X[n1]; mn[1] += X[n1+1]; n1+=2; }
        cblas_zaxpy(N,mn,ni,0,X,1);
        cblas_zdotc_sub(N,X,1,X,1,(_Complex double *)sm2);
        for (n=0; n<2*N; n+=2)
        {
            x1[n] = X[n]*X[n] - X[n+1]*X[n+1];
            x1[n+1] = X[n]*X[n+1] + X[n+1]*X[n];
        }
        for (n=0; n<2*N; n+=2)
        {
            tmp = x1[n]*X[n] - x1[n+1]*X[n+1];
            x1[n+1] = x1[n]*X[n+1] + x1[n+1]*X[n];
            x1[n] = tmp;
        }
        cblas_zdotu_sub(N,x1,1,X,1,(_Complex double *)sm3);
        Y[0] = N*sm3[0]/(sm2[0]*sm2[0]); Y[1] = N*sm3[1]/(sm2[0]*sm2[0]);
        if (!biased)
        {
            Y[0] =  3.0 + (Y[0]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3));
            Y[1] *= (N1+1)*(N1-1)/((N1-2)*(N1-3));
        }
    }
    else
    {
        if (iscolmajor) {}
        fprintf(stderr,"error in skewness_z: complex skewness only supported for vectors currently \n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
