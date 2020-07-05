//Gets kurtosis of each row or col of X according to dim.
//For complex case, output is complex.
//I follow the Octave convention for complex kurtosis (but see literature for other ideas later).
//Complex kurtosis only supported for vectors currently.
//This works in place.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>
//struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
//clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int kurtosis_s (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased);
int kurtosis_d (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased);
int kurtosis_c (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased);
int kurtosis_z (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased);


int kurtosis_s (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    if (R<1) { fprintf(stderr,"error in kurtosis_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in kurtosis_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in kurtosis_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in kurtosis_s: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int N2 = N/N1;
    const float o = 1.0f, ni = 1.0f / N1;

    if (N1<2) { fprintf(stderr,"error in kurtosis_s: N must be > 1 \n"); return 1; }
    else if (N1==N)
    {
        float sm2 = 0.0f, sm4 = 0.0f;
        Y[0] = 0.0f;
        for (int n1=0; n1<N1; n1++) { Y[0] += X[n1]; }
        cblas_saxpy(N1,-Y[0],&ni,0,X,1);
        for (int n=0; n<N; n++) { X[n] *= X[n]; sm2 += X[n]; sm4 += X[n]*X[n]; }
        Y[0] = N * sm4 / (sm2*sm2);
        if (!biased) { Y[0] =  3.0f + (Y[0]*(N+1)-3*(N-1)) * (N-1)/((N-2)*(N-3)); }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : N1;
        const int ldb = (iscolmajor) ? RC/N1 : 1;
        const int ldc = (iscolmajor) ? R : C;
        float *x1;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N1,&ni,0,x1,1);
        cblas_sgemv(Ord,Tr,(int)R,(int)C,1.0f,X,ldc,x1,1,0.0f,Y,1);
        cblas_scopy(N1,&o,0,x1,1);
        if (dim==0) { cblas_sgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,x1,lda,Y,ldb,1.0f,X,ldc); }
        else { cblas_sgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,Y,ldb,x1,lda,1.0f,X,ldc); }
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (int n=0; n<N; n++) { X[n] *= X[n]; }
            for (int n2=0; n2<N2; n2++)
            {
                Y[n2] = cblas_sdot(N1,&X[n2*N1],1,x1,1);
                Y[n2] = N1 * cblas_sdot(N1,&X[n2*N1],1,&X[n2*N1],1) / (Y[n2]*Y[n2]);
            }
        }
        else
        {
            float *sm2, *sm4;
            if (!(sm2=(float *)calloc((size_t)N2,sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm4=(float *)calloc((size_t)N2,sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with calloc. "); perror("calloc"); return 1; }
            for (int n=0; n<N; n++) { int n2 = n%N2; X[n] *= X[n]; sm2[n2] += X[n]; sm4[n2] += X[n]*X[n]; }
            for (int n2=0; n2<N2; n2++) { Y[n2] = N1 * sm4[n2] / (sm2[n2]*sm2[n2]); }
            free(sm2); free(sm4);
        }
        if (!biased)
        {
            for (int n2=0; n2<N2; n2++) { Y[n2] =  3.0f + (Y[n2]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        float *x1, *xni;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(N1,&o,0,x1,1); cblas_scopy(N1,&ni,0,xni,1);
        for (size_t r=0, n=0, n2=0; r<R; r++, n+=C*S, n2+=C)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0f,&X[n],(int)S,xni,1,0.0f,&Y[n2],1);
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)C,(int)S,1,-1.0f,&Y[n2],1,x1,(int)S,1.0f,&X[n],(int)S);
        }
        for (int n=0; n<N; n++) { X[n] *= X[n]; }
        for (size_t r=0, n=0, n2=0; r<R; r++, n+=C*S, n2+=C)
        {
            for (size_t c=0; c<C; c++)
            {
                Y[n2+c] = cblas_sdot(N1,&X[n+c*S],1,x1,1);
                Y[n2+c] = N1 * cblas_sdot(N1,&X[n+c*S],1,&X[n+c*S],1) / (Y[n2+c]*Y[n2+c]);
            }
        }
        if (!biased)
        {
            for (int n2=0; n2<N2; n2++) { Y[n2] =  3.0f + (Y[n2]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1); free(xni);
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        float sm1, sm2, sm4, *x1;
        if (!(x1=(float *)malloc((size_t)N1*sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with malloc. "); perror("malloc"); return 1; }
        for (int l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
        {
            for (int m=0; m<M; m++, n+=J, n2++)
            {
                cblas_scopy(N1,&X[n],K,x1,1);
                sm1 = cblas_sdot(N1,x1,1,&o,0);
                cblas_saxpy(N1,-sm1,&ni,0,x1,1);
                sm2 = sm4 = 0.0f;
                for (int n1=0; n1<N1; n1++) { float x = x1[n1]*x1[n1]; sm2 += x; sm4 += x*x; }
                Y[n2] = N1 * sm4 / (sm2*sm2);
            }
        }
        if (!biased)
        {
            for (int n=0; n<N/N1; n++) { Y[n] = 3.0f + (Y[n]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1);
    }

    return 0;
}


int kurtosis_d (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    if (R<1) { fprintf(stderr,"error in kurtosis_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in kurtosis_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in kurtosis_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in kurtosis_d: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int N2 = N/N1;
    const double o = 1.0, ni = 1.0 / N1;

    if (N1<2) { fprintf(stderr,"error in kurtosis_d: N must be > 1 \n"); return 1; }
    else if (N1==N)
    {
        double sm2 = 0.0, sm4 = 0.0;
        Y[0] = 0.0;
        for (int n1=0; n1<N1; n1++) { Y[0] += X[n1]; }
        cblas_daxpy(N1,-Y[0],&ni,0,X,1);
        for (int n=0; n<N; n++) { X[n] *= X[n]; sm2 += X[n]; sm4 += X[n]*X[n]; }
        Y[0] = N * sm4 / (sm2*sm2);
        if (!biased) { Y[0] =  3.0 + (Y[0]*(N+1)-3*(N-1)) * (N-1)/((N-2)*(N-3)); }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : N1;
        const int ldb = (iscolmajor) ? RC/N1 : 1;
        const int ldc = (iscolmajor) ? R : C;
        double *x1;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N1,&ni,0,x1,1);
        cblas_dgemv(Ord,Tr,(int)R,(int)C,1.0,X,ldc,x1,1,0.0,Y,1);
        cblas_dcopy(N1,&o,0,x1,1);
        if (dim==0) { cblas_dgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0,x1,lda,Y,ldb,1.0,X,ldc); }
        else { cblas_dgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0,Y,ldb,x1,lda,1.0,X,ldc); }
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (int n=0; n<N; n++) { X[n] *= X[n]; }
            for (int n2=0; n2<N2; n2++)
            {
                Y[n2] = cblas_ddot(N1,&X[n2*N1],1,x1,1);
                Y[n2] = N1 * cblas_ddot(N1,&X[n2*N1],1,&X[n2*N1],1) / (Y[n2]*Y[n2]);
            }
        }
        else
        {
            double *sm2, *sm4;
            if (!(sm2=(double *)calloc((size_t)N2,sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm4=(double *)calloc((size_t)N2,sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with calloc. "); perror("calloc"); return 1; }
            for (int n=0; n<N; n++) { int n2 = n%N2; X[n] *= X[n]; sm2[n2] += X[n]; sm4[n2] += X[n]*X[n]; }
            for (int n2=0; n2<N2; n2++) { Y[n2] = N1 * sm4[n2] / (sm2[n2]*sm2[n2]); }
            free(sm2); free(sm4);
        }
        if (!biased)
        {
            for (int n2=0; n2<N2; n2++) { Y[n2] =  3.0 + (Y[n2]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        double *x1, *xni;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(N1,&o,0,x1,1); cblas_dcopy(N1,&ni,0,xni,1);
        for (size_t r=0, n=0, n2=0; r<R; r++, n+=C*S, n2+=C)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0,&X[n],(int)S,xni,1,0.0,&Y[n2],1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)C,(int)S,1,-1.0,&Y[n2],1,x1,(int)S,1.0,&X[n],(int)S);
        }
        for (int n=0; n<N; n++) { X[n] *= X[n]; }
        for (size_t r=0, n=0, n2=0; r<R; r++, n+=C*S, n2+=C)
        {
            for (size_t c=0; c<C; c++)
            {
                Y[n2+c] = cblas_ddot(N1,&X[n+c*S],1,x1,1);
                Y[n2+c] = N1 * cblas_ddot(N1,&X[n+c*S],1,&X[n+c*S],1) / (Y[n2+c]*Y[n2+c]);
            }
        }
        if (!biased)
        {
            for (int n2=0; n2<N2; n2++) { Y[n2] = 3.0 + (Y[n2]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1); free(xni);
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        double sm1, sm2, sm4, *x1;
        if (!(x1=(double *)malloc((size_t)N1*sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with malloc. "); perror("malloc"); return 1; }
        for (int l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
        {
            for (int m=0; m<M; m++, n+=J, n2++)
            {
                cblas_dcopy(N1,&X[n],K,x1,1);
                sm1 = cblas_ddot(N1,x1,1,&o,0);
                cblas_daxpy(N1,-sm1,&ni,0,x1,1);
                sm2 = sm4 = 0.0;
                for (int n1=0; n1<N1; n1++) { double x = x1[n1]*x1[n1]; sm2 += x; sm4 += x*x; }
                Y[n2] = N1 * sm4 / (sm2*sm2);
            }
        }
        if (!biased)
        {
            for (int n=0; n<N/N1; n++) { Y[n] =  3.0 + (Y[n]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1);
    }

    return 0;
}


int kurtosis_c (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    if (R<1) { fprintf(stderr,"error in kurtosis_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in kurtosis_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in kurtosis_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in kurtosis_c: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ni[2] = {-1.0f/N1,0.0f};

    if (N1==1)
    {
        const float z[2] =  {0.0f,0.0f};
        cblas_ccopy(N,z,0,Y,1);
    }
    else if (N1==N)
    {
        float sm1[2] = {0.0f,0.0f}, sm2[2] = {0.0f,0.0f}, sm4[2] = {0.0f,0.0f};
        float tmp, *x2;
        if (!(x2=(float *)malloc((size_t)(2*N1)*sizeof(float)))) { fprintf(stderr,"error in kurtosis_c: problem with malloc. "); perror("malloc"); return 1; }
        for (int n=0; n<2*N; n+=2) { sm1[0] += X[n]; sm1[1] += X[n+1]; }
        cblas_caxpy(N,sm1,ni,0,X,1);
        cblas_cdotc_sub(N,X,1,X,1,(_Complex float *)sm2);
        for (int n=0; n<2*N; n+=2)
        {
            x2[n] = X[n]*X[n] - X[n+1]*X[n+1];
            x2[n+1] = X[n]*X[n+1] + X[n+1]*X[n];
            tmp = x2[n]*X[n] - x2[n+1]*X[n+1];
            x2[n+1] = x2[n]*X[n+1] + x2[n+1]*X[n];
            x2[n] = tmp;
        }
        cblas_cdotu_sub(N,x2,1,X,1,(_Complex float *)sm4);
        Y[0] = N*sm4[0]/(sm2[0]*sm2[0]); Y[1] = N*sm4[1]/(sm2[0]*sm2[0]);
        if (!biased)
        {
            Y[0] = 3.0f + (Y[0]*(N+1)-3*(N-1)) * (N-1)/((N-2)*(N-3));
            Y[1] *= (N+1)*(N-1) / (float)((N-2)*(N-3));
        }
        free(x2);
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        float tmp, sm1[2], sm2[2], sm4[2], *x1, *x2;
        if (!(x1=(float *)malloc((size_t)(2*N1)*sizeof(float)))) { fprintf(stderr,"error in kurtosis_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(x2=(float *)malloc((size_t)(2*N1)*sizeof(float)))) { fprintf(stderr,"error in kurtosis_c: problem with malloc. "); perror("malloc"); return 1; }
        for (int l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
        {
            for (int m=0; m<M; m++, n+=2*J, n2+=2)
            {
                cblas_ccopy(N1,&X[n],K,x1,1);
                sm1[0] = sm1[1] = sm2[0] = sm2[1] = sm4[0] = sm4[1] = 0.0f;
                for (int n1=0; n1<2*N1; n1+=2) { sm1[0] += x1[n1]; sm1[1] += x1[n1+1]; }
                cblas_caxpy(N1,sm1,ni,0,x1,1);
                cblas_cdotc_sub(N1,x1,1,x1,1,(_Complex float *)sm2);
                for (int n1=0; n1<2*N1; n1+=2)
                {
                    x2[n1] = x1[n1]*x1[n1] - x1[n1+1]*x1[n1+1];
                    x2[n1+1] = x1[n1]*x1[n1+1] + x1[n1+1]*x1[n1];
                    tmp = x2[n1]*x1[n1] - x2[n1+1]*x1[n1+1];
                    x2[n1+1] = x2[n1]*x1[n1+1] + x2[n1+1]*x1[n1];
                    x2[n1] = tmp;
                }
                cblas_cdotu_sub(N1,x2,1,x1,1,(_Complex float *)sm4);
                Y[n2] = N1*sm4[0]/(sm2[0]*sm2[0]); Y[n2+1] = N1*sm4[1]/(sm2[0]*sm2[0]);
            }
        }
        if (!biased)
        {
            for (int n2=0; n2<2*N/N1; n2+=2)
            {
                Y[n2] = 3.0f + (Y[n2]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3));
                Y[n2+1] *= (N1+1)*(N1-1) / (float)((N1-2)*(N1-3));
            }
        }
        free(x1); free(x2);
    }

    return 0;
}


int kurtosis_z (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    if (R<1) { fprintf(stderr,"error in kurtosis_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in kurtosis_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in kurtosis_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in kurtosis_z: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni[2] = {-1.0/N1,0.0};

    if (N1==1)
    {
        const double z[2] =  {0.0,0.0};
        cblas_zcopy(N,z,0,Y,1);
    }
    else if (N1==N)
    {
        double sm1[2] = {0.0,0.0}, sm2[2] = {0.0,0.0}, sm4[2] = {0.0,0.0};
        double tmp, *x2;
        if (!(x2=(double *)malloc((size_t)(2*N1)*sizeof(double)))) { fprintf(stderr,"error in kurtosis_z: problem with malloc. "); perror("malloc"); return 1; }
        for (int n=0; n<2*N; n+=2) { sm1[0] += X[n]; sm1[1] += X[n+1]; }
        cblas_zaxpy(N,sm1,ni,0,X,1);
        cblas_zdotc_sub(N,X,1,X,1,(_Complex double *)sm2);
        for (int n=0; n<2*N; n+=2)
        {
            x2[n] = X[n]*X[n] - X[n+1]*X[n+1];
            x2[n+1] = X[n]*X[n+1] + X[n+1]*X[n];
            tmp = x2[n]*X[n] - x2[n+1]*X[n+1];
            x2[n+1] = x2[n]*X[n+1] + x2[n+1]*X[n];
            x2[n] = tmp;
        }
        cblas_zdotu_sub(N,x2,1,X,1,(_Complex double *)sm4);
        Y[0] = N*sm4[0]/(sm2[0]*sm2[0]); Y[1] = N*sm4[1]/(sm2[0]*sm2[0]);
        if (!biased)
        {
            Y[0] =  3.0 + (Y[0]*(N+1)-3*(N-1)) * (N-1)/((N-2)*(N-3));
            Y[1] *= (N+1)*(N-1) / (double)((N-2)*(N-3));
        }
        free(x2);
    }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int L = N/(M*N1);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        double tmp, sm1[2], sm2[2], sm4[2], *x1, *x2;
        if (!(x1=(double *)malloc((size_t)(2*N1)*sizeof(double)))) { fprintf(stderr,"error in kurtosis_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(x2=(double *)malloc((size_t)(2*N1)*sizeof(double)))) { fprintf(stderr,"error in kurtosis_z: problem with malloc. "); perror("malloc"); return 1; }
        for (int l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
        {
            for (int m=0; m<M; m++, n+=2*J, n2+=2)
            {
                cblas_zcopy(N1,&X[n],K,x1,1);
                sm1[0] = sm1[1] = sm2[0] = sm2[1] = sm4[0] = sm4[1] = 0.0;
                for (int n1=0; n1<2*N1; n1+=2) { sm1[0] += x1[n1]; sm1[1] += x1[n1+1]; }
                cblas_zaxpy(N1,sm1,ni,0,x1,1);
                cblas_zdotc_sub(N1,x1,1,x1,1,(_Complex double *)sm2);
                for (int n1=0; n1<2*N1; n1+=2)
                {
                    x2[n1] = x1[n1]*x1[n1] - x1[n1+1]*x1[n1+1];
                    x2[n1+1] = x1[n1]*x1[n1+1] + x1[n1+1]*x1[n1];
                    tmp = x2[n1]*x1[n1] - x2[n1+1]*x1[n1+1];
                    x2[n1+1] = x2[n1]*x1[n1+1] + x2[n1+1]*x1[n1];
                    x2[n1] = tmp;
                }
                cblas_zdotu_sub(N1,x2,1,x1,1,(_Complex double *)sm4);
                Y[n2] = N1*sm4[0]/(sm2[0]*sm2[0]); Y[n2+1] = N1*sm4[1]/(sm2[0]*sm2[0]);
            }
        }
        if (!biased)
        {
            for (int n2=0; n2<2*N/N1; n2+=2)
            {
                Y[n2] = 3.0 + (Y[n2]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3));
                Y[n2+1] *= (N1+1)*(N1-1) / (double)((N1-2)*(N1-3));
            }
        }
        free(x1); free(x2);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
