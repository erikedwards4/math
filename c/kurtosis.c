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
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t N2 = N/N1;
    const float o = 1.0f, ni = 1.0f / N1;

    if (N1<4) { fprintf(stderr,"error in kurtosis_s: N1 must be > 3 \n"); return 1; }
    else if (N1==N)
    {
        float sm2 = 0.0f, sm4 = 0.0f;
        Y[0] = 0.0f;
        for (size_t n1=0; n1<N1; n1++) { Y[0] += X[n1]; }
        cblas_saxpy((int)N1,-Y[0],&ni,0,X,1);
        for (size_t n=0; n<N; n++) { X[n] *= X[n]; sm2 += X[n]; sm4 += X[n]*X[n]; }
        Y[0] = N * sm4 / (sm2*sm2);
        if (!biased) { Y[0] =  3.0f + (Y[0]*(N+1)-3*(N-1)) * (N-1)/((N-2)*(N-3)); }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : (int)N1;
        const int ldb = (iscolmajor) ? (int)N2 : 1;
        const int ldc = (iscolmajor) ? (int)R : (int)C;
        float *x1;
        if (!(x1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)N1,&ni,0,x1,1);
        cblas_sgemv(Ord,Tr,(int)R,(int)C,1.0f,X,ldc,x1,1,0.0f,Y,1);
        cblas_scopy((int)N1,&o,0,x1,1);
        if (dim==0) { cblas_sgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,x1,lda,Y,ldb,1.0f,X,ldc); }
        else { cblas_sgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,Y,ldb,x1,lda,1.0f,X,ldc); }
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
            for (size_t n2=0; n2<N2; n2++, X+=N1)
            {
                Y[n2] = cblas_sdot((int)N1,X,1,x1,1);
                Y[n2] = N1 * cblas_sdot((int)N1,X,1,X,1) / (Y[n2]*Y[n2]);
            }
        }
        else
        {
            float *sm2, *sm4;
            if (!(sm2=(float *)calloc((size_t)N2,sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm4=(float *)calloc((size_t)N2,sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t n=0; n<N; n++) { size_t n2 = n%N2; X[n] *= X[n]; sm2[n2] += X[n]; sm4[n2] += X[n]*X[n]; }
            for (size_t n2=0; n2<N2; n2++) { Y[n2] = N1 * sm4[n2] / (sm2[n2]*sm2[n2]); }
            free(sm2); free(sm4);
        }
        if (!biased)
        {
            for (size_t n2=0; n2<N2; n2++) { Y[n2] =  3.0f + (Y[n2]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        float *x1, *xni;
        if (!(x1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)N1,&o,0,x1,1); cblas_scopy((int)N1,&ni,0,xni,1);
        for (size_t r=0; r<R; r++, X+=C*S, Y+=C)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0f,X,(int)S,xni,1,0.0f,Y,1);
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)C,(int)S,1,-1.0f,Y,1,x1,(int)S,1.0f,X,(int)S);
        }
        X -= N; Y -= RC;
        for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, X+=S, Y++)
            {
                *Y = cblas_sdot((int)N1,X,1,x1,1);
                *Y = N1 * cblas_sdot((int)N1,X,1,X,1) / (*Y**Y);
            }
        }
        if (!biased)
        {
            Y -= RC;
            for (size_t n2=0; n2<N2; n2++) { Y[n2] =  3.0f + (Y[n2]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1); free(xni);
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        float sm1, sm2, sm4, *x1;
        if (!(x1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J, Y++)
            {
                cblas_scopy((int)N1,X,(int)K,x1,1);
                sm1 = cblas_sdot((int)N1,x1,1,&o,0);
                cblas_saxpy((int)N1,-sm1,&ni,0,x1,1);
                sm2 = sm4 = 0.0f;
                for (size_t n1=0; n1<N1; n1++) { float x = x1[n1]*x1[n1]; sm2 += x; sm4 += x*x; }
                *Y = N1 * sm4 / (sm2*sm2);
            }
        }
        if (!biased)
        {
            Y -= L*M;
            for (size_t n=0; n<N/N1; n++) { Y[n] = 3.0f + (Y[n]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1);
    }

    return 0;
}


int kurtosis_d (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t N2 = N/N1;
    const double o = 1.0, ni = 1.0 / N1;

    if (N1<4) { fprintf(stderr,"error in kurtosis_d: N1 must be > 3 \n"); return 1; }
    else if (N1==N)
    {
        double sm2 = 0.0, sm4 = 0.0;
        Y[0] = 0.0;
        for (size_t n1=0; n1<N1; n1++) { Y[0] += X[n1]; }
        cblas_daxpy((int)N1,-Y[0],&ni,0,X,1);
        for (size_t n=0; n<N; n++) { X[n] *= X[n]; sm2 += X[n]; sm4 += X[n]*X[n]; }
        Y[0] = N * sm4 / (sm2*sm2);
        if (!biased) { Y[0] =  3.0 + (Y[0]*(N+1)-3*(N-1)) * (N-1)/((N-2)*(N-3)); }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : (int)N1;
        const int ldb = (iscolmajor) ? (int)N2 : 1;
        const int ldc = (iscolmajor) ? (int)R : (int)C;
        double *x1;
        if (!(x1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)N1,&ni,0,x1,1);
        cblas_dgemv(Ord,Tr,(int)R,(int)C,1.0,X,ldc,x1,1,0.0,Y,1);
        cblas_dcopy((int)N1,&o,0,x1,1);
        if (dim==0) { cblas_dgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0,x1,lda,Y,ldb,1.0,X,ldc); }
        else { cblas_dgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0,Y,ldb,x1,lda,1.0,X,ldc); }
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
            for (size_t n2=0; n2<N2; n2++, X+=N1)
            {
                Y[n2] = cblas_ddot((int)N1,X,1,x1,1);
                Y[n2] = N1 * cblas_ddot((int)N1,X,1,X,1) / (Y[n2]*Y[n2]);
            }
        }
        else
        {
            double *sm2, *sm4;
            if (!(sm2=(double *)calloc((size_t)N2,sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm4=(double *)calloc((size_t)N2,sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t n=0; n<N; n++) { size_t n2 = n%N2; X[n] *= X[n]; sm2[n2] += X[n]; sm4[n2] += X[n]*X[n]; }
            for (size_t n2=0; n2<N2; n2++) { Y[n2] = N1 * sm4[n2] / (sm2[n2]*sm2[n2]); }
            free(sm2); free(sm4);
        }
        if (!biased)
        {
            for (size_t n2=0; n2<N2; n2++) { Y[n2] =  3.0 + (Y[n2]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        double *x1, *xni;
        if (!(x1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)N1,&o,0,x1,1); cblas_dcopy((int)N1,&ni,0,xni,1);
        for (size_t r=0; r<R; r++, X+=C*S, Y+=C)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)C,(int)S,1.0,X,(int)S,xni,1,0.0,Y,1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)C,(int)S,1,-1.0,Y,1,x1,(int)S,1.0,X,(int)S);
        }
        X -= N; Y -= RC;
        for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
        for (size_t r=0; r<R; r++)
        {
            for (size_t c=0; c<C; c++, X+=S, Y++)
            {
                *Y = cblas_ddot((int)N1,X,1,x1,1);
                *Y = N1 * cblas_ddot((int)N1,X,1,X,1) / (*Y**Y);
            }
        }
        if (!biased)
        {
            Y -= RC;
            for (size_t n2=0; n2<N2; n2++) { Y[n2] =  3.0 + (Y[n2]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1); free(xni);
    }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        double sm1, sm2, sm4, *x1;
        if (!(x1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J, Y++)
            {
                cblas_dcopy((int)N1,X,(int)K,x1,1);
                sm1 = cblas_ddot((int)N1,x1,1,&o,0);
                cblas_daxpy((int)N1,-sm1,&ni,0,x1,1);
                sm2 = sm4 = 0.0;
                for (size_t n1=0; n1<N1; n1++) { double x = x1[n1]*x1[n1]; sm2 += x; sm4 += x*x; }
                *Y = N1 * sm4 / (sm2*sm2);
            }
        }
        if (!biased)
        {
            Y -= L*M;
            for (size_t n=0; n<N/N1; n++) { Y[n] =  3.0 + (Y[n]*(N1+1)-3*(N1-1)) * (N1-1)/((N1-2)*(N1-3)); }
        }
        free(x1);
    }

    return 0;
}


int kurtosis_c (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ni[2] = {-1.0f/N1,0.0f};

    if (N1<4)  { fprintf(stderr,"error in kurtosis_c: N1 must be > 3 \n"); return 1; }
    else if (N1==N)
    {
        float sm1[2] = {0.0f,0.0f}, sm2[2] = {0.0f,0.0f}, sm4[2] = {0.0f,0.0f};
        float tmp, *x2;
        if (!(x2=(float *)malloc(2*N1*sizeof(float)))) { fprintf(stderr,"error in kurtosis_c: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t n=0; n<2*N; n+=2) { sm1[0] += X[n]; sm1[1] += X[n+1]; }
        cblas_caxpy((int)N,sm1,ni,0,X,1);
        cblas_cdotc_sub((int)N,X,1,X,1,(_Complex float *)sm2);
        for (size_t n=0; n<2*N; n+=2)
        {
            x2[n] = X[n]*X[n] - X[n+1]*X[n+1];
            x2[n+1] = X[n]*X[n+1] + X[n+1]*X[n];
            tmp = x2[n]*X[n] - x2[n+1]*X[n+1];
            x2[n+1] = x2[n]*X[n+1] + x2[n+1]*X[n];
            x2[n] = tmp;
        }
        cblas_cdotu_sub((int)N,x2,1,X,1,(_Complex float *)sm4);
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
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        float tmp, sm1[2], sm2[2], sm4[2], *x1, *x2;
        if (!(x1=(float *)malloc(2*N1*sizeof(float)))) { fprintf(stderr,"error in kurtosis_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(x2=(float *)malloc(2*N1*sizeof(float)))) { fprintf(stderr,"error in kurtosis_c: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++, X+=2*M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=2*J)
            {
                cblas_ccopy((int)N1,X,(int)K,x1,1);
                sm1[0] = sm1[1] = sm2[0] = sm2[1] = sm4[0] = sm4[1] = 0.0f;
                for (size_t n1=0; n1<2*N1; n1+=2) { sm1[0] += x1[n1]; sm1[1] += x1[n1+1]; }
                cblas_caxpy((int)N1,sm1,ni,0,x1,1);
                cblas_cdotc_sub((int)N1,x1,1,x1,1,(_Complex float *)sm2);
                for (size_t n1=0; n1<2*N1; n1+=2)
                {
                    x2[n1] = x1[n1]*x1[n1] - x1[n1+1]*x1[n1+1];
                    x2[n1+1] = x1[n1]*x1[n1+1] + x1[n1+1]*x1[n1];
                    tmp = x2[n1]*x1[n1] - x2[n1+1]*x1[n1+1];
                    x2[n1+1] = x2[n1]*x1[n1+1] + x2[n1+1]*x1[n1];
                    x2[n1] = tmp;
                }
                cblas_cdotu_sub((int)N1,x2,1,x1,1,(_Complex float *)sm4);
                *Y++ = N1*sm4[0]/(sm2[0]*sm2[0]);
                *Y++ = N1*sm4[1]/(sm2[0]*sm2[0]);
            }
        }
        if (!biased)
        {
            Y -= 2*L*M;
            for (size_t n2=0; n2<2*N/N1; n2+=2)
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
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni[2] = {-1.0/N1,0.0};

    if (N1<4)  { fprintf(stderr,"error in kurtosis_z: N1 must be > 3 \n"); return 1; }
    else if (N1==N)
    {
        double sm1[2] = {0.0,0.0}, sm2[2] = {0.0,0.0}, sm4[2] = {0.0,0.0};
        double tmp, *x2;
        if (!(x2=(double *)malloc(2*N1*sizeof(double)))) { fprintf(stderr,"error in kurtosis_z: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t n=0; n<2*N; n+=2) { sm1[0] += X[n]; sm1[1] += X[n+1]; }
        cblas_zaxpy((int)N,sm1,ni,0,X,1);
        cblas_zdotc_sub((int)N,X,1,X,1,(_Complex double *)sm2);
        for (size_t n=0; n<2*N; n+=2)
        {
            x2[n] = X[n]*X[n] - X[n+1]*X[n+1];
            x2[n+1] = X[n]*X[n+1] + X[n+1]*X[n];
            tmp = x2[n]*X[n] - x2[n+1]*X[n+1];
            x2[n+1] = x2[n]*X[n+1] + x2[n+1]*X[n];
            x2[n] = tmp;
        }
        cblas_zdotu_sub((int)N,x2,1,X,1,(_Complex double *)sm4);
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
        const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t L = N/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        double tmp, sm1[2], sm2[2], sm4[2], *x1, *x2;
        if (!(x1=(double *)malloc(2*N1*sizeof(double)))) { fprintf(stderr,"error in kurtosis_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(x2=(double *)malloc(2*N1*sizeof(double)))) { fprintf(stderr,"error in kurtosis_z: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++, X+=2*M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=2*J)
            {
                cblas_zcopy((int)N1,X,(int)K,x1,1);
                sm1[0] = sm1[1] = sm2[0] = sm2[1] = sm4[0] = sm4[1] = 0.0;
                for (size_t n1=0; n1<2*N1; n1+=2) { sm1[0] += x1[n1]; sm1[1] += x1[n1+1]; }
                cblas_zaxpy((int)N1,sm1,ni,0,x1,1);
                cblas_zdotc_sub((int)N1,x1,1,x1,1,(_Complex double *)sm2);
                for (size_t n1=0; n1<2*N1; n1+=2)
                {
                    x2[n1] = x1[n1]*x1[n1] - x1[n1+1]*x1[n1+1];
                    x2[n1+1] = x1[n1]*x1[n1+1] + x1[n1+1]*x1[n1];
                    tmp = x2[n1]*x1[n1] - x2[n1+1]*x1[n1+1];
                    x2[n1+1] = x2[n1]*x1[n1+1] + x2[n1+1]*x1[n1];
                    x2[n1] = tmp;
                }
                cblas_zdotu_sub((int)N1,x2,1,x1,1,(_Complex double *)sm4);
                *Y++ = N1*sm4[0]/(sm2[0]*sm2[0]);
                *Y++ = N1*sm4[1]/(sm2[0]*sm2[0]);
            }
        }
        if (!biased)
        {
            Y -= 2*L*M;
            for (size_t n2=0; n2<2*N/N1; n2+=2)
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
