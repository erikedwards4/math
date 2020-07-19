//Gets kurtosis of each row or col of X according to dim.
//For complex case, output is complex.
//I follow the Octave convention for complex kurtosis (but see literature for other ideas later).
//Complex kurtosis only supported for vectors currently.
//This works in place.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int kurtosis_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int kurtosis_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int kurtosis_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int kurtosis_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);


int kurtosis_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in kurtosis_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t V = N/L;
    const float o = 1.0f, ni = 1.0f / L;

    if (L<4) { fprintf(stderr,"error in kurtosis_s: L must be > 3 \n"); return 1; }
    else if (L==N)
    {
        float sm2 = 0.0f, sm4 = 0.0f;
        Y[0] = 0.0f;
        for (size_t l=0; l<L; l++) { Y[0] += X[l]; }
        cblas_saxpy((int)L,-Y[0],&ni,0,X,1);
        for (size_t n=0; n<N; n++) { X[n] *= X[n]; sm2 += X[n]; sm4 += X[n]*X[n]; }
        Y[0] = N * sm4 / (sm2*sm2);
        if (!biased) { Y[0] =  3.0f + (Y[0]*(N+1)-3*(N-1)) * (N-1)/((N-2)*(N-3)); }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : (int)L;
        const int ldb = (iscolmajor) ? (int)V : 1;
        const int ldc = (iscolmajor) ? (int)R : (int)C;
        float *x1;
        if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)L,&ni,0,x1,1);
        cblas_sgemv(Ord,Tr,(int)R,(int)C,1.0f,X,ldc,x1,1,0.0f,Y,1);
        cblas_scopy((int)L,&o,0,x1,1);
        if (dim==0) { cblas_sgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,x1,lda,Y,ldb,1.0f,X,ldc); }
        else { cblas_sgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0f,Y,ldb,x1,lda,1.0f,X,ldc); }
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
            for (size_t v=0; v<V; v++, X+=L)
            {
                Y[v] = cblas_sdot((int)L,X,1,x1,1);
                Y[v] = L * cblas_sdot((int)L,X,1,X,1) / (Y[v]*Y[v]);
            }
        }
        else
        {
            float *sm2, *sm4;
            if (!(sm2=(float *)calloc((size_t)V,sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm4=(float *)calloc((size_t)V,sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t n=0; n<N; n++) { size_t v = n%V; X[n] *= X[n]; sm2[v] += X[n]; sm4[v] += X[n]*X[n]; }
            for (size_t v=0; v<V; v++) { Y[v] = L * sm4[v] / (sm2[v]*sm2[v]); }
            free(sm2); free(sm4);
        }
        if (!biased)
        {
            for (size_t v=0; v<V; v++) { Y[v] =  3.0f + (Y[v]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        float *x1, *xni;
        if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)L,&o,0,x1,1); cblas_scopy((int)L,&ni,0,xni,1);
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
                *Y = cblas_sdot((int)L,X,1,x1,1);
                *Y = L * cblas_sdot((int)L,X,1,X,1) / (*Y**Y);
            }
        }
        if (!biased)
        {
            Y -= RC;
            for (size_t v=0; v<V; v++) { Y[v] =  3.0f + (Y[v]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
        }
        free(x1); free(xni);
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        float sm1, sm2, sm4, *x1;
        if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=J, Y++)
            {
                cblas_scopy((int)L,X,(int)K,x1,1);
                sm1 = cblas_sdot((int)L,x1,1,&o,0);
                cblas_saxpy((int)L,-sm1,&ni,0,x1,1);
                sm2 = sm4 = 0.0f;
                for (size_t l=0; l<L; l++) { float x = x1[l]*x1[l]; sm2 += x; sm4 += x*x; }
                *Y = L * sm4 / (sm2*sm2);
            }
        }
        if (!biased)
        {
            Y -= G*B;
            for (size_t n=0; n<N/L; n++) { Y[n] = 3.0f + (Y[n]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
        }
        free(x1);
    }

    return 0;
}


int kurtosis_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in kurtosis_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t V = N/L;
    const double o = 1.0, ni = 1.0 / L;

    if (L<4) { fprintf(stderr,"error in kurtosis_d: L must be > 3 \n"); return 1; }
    else if (L==N)
    {
        double sm2 = 0.0, sm4 = 0.0;
        Y[0] = 0.0;
        for (size_t l=0; l<L; l++) { Y[0] += X[l]; }
        cblas_daxpy((int)L,-Y[0],&ni,0,X,1);
        for (size_t n=0; n<N; n++) { X[n] *= X[n]; sm2 += X[n]; sm4 += X[n]*X[n]; }
        Y[0] = N * sm4 / (sm2*sm2);
        if (!biased) { Y[0] =  3.0 + (Y[0]*(N+1)-3*(N-1)) * (N-1)/((N-2)*(N-3)); }
    }
    else if (SH==1)
    {
        const CBLAS_TRANSPOSE Tr = (dim==0) ? CblasTrans : CblasNoTrans;
        const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
        const int lda = (iscolmajor) ? 1 : (int)L;
        const int ldb = (iscolmajor) ? (int)V : 1;
        const int ldc = (iscolmajor) ? (int)R : (int)C;
        double *x1;
        if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)L,&ni,0,x1,1);
        cblas_dgemv(Ord,Tr,(int)R,(int)C,1.0,X,ldc,x1,1,0.0,Y,1);
        cblas_dcopy((int)L,&o,0,x1,1);
        if (dim==0) { cblas_dgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0,x1,lda,Y,ldb,1.0,X,ldc); }
        else { cblas_dgemm(Ord,Tr,Tr,(int)R,(int)C,1,-1.0,Y,ldb,x1,lda,1.0,X,ldc); }
        if ((dim==0 && iscolmajor) || (dim==1 && !iscolmajor))
        {
            for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
            for (size_t v=0; v<V; v++, X+=L)
            {
                Y[v] = cblas_ddot((int)L,X,1,x1,1);
                Y[v] = L * cblas_ddot((int)L,X,1,X,1) / (Y[v]*Y[v]);
            }
        }
        else
        {
            double *sm2, *sm4;
            if (!(sm2=(double *)calloc((size_t)V,sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm4=(double *)calloc((size_t)V,sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t n=0; n<N; n++) { size_t v = n%V; X[n] *= X[n]; sm2[v] += X[n]; sm4[v] += X[n]*X[n]; }
            for (size_t v=0; v<V; v++) { Y[v] = L * sm4[v] / (sm2[v]*sm2[v]); }
            free(sm2); free(sm4);
        }
        if (!biased)
        {
            for (size_t v=0; v<V; v++) { Y[v] =  3.0 + (Y[v]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
        }
        free(x1);
    }
    else if (!iscolmajor && dim==2 && H==1)
    {
        double *x1, *xni;
        if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(xni=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)L,&o,0,x1,1); cblas_dcopy((int)L,&ni,0,xni,1);
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
                *Y = cblas_ddot((int)L,X,1,x1,1);
                *Y = L * cblas_ddot((int)L,X,1,X,1) / (*Y**Y);
            }
        }
        if (!biased)
        {
            Y -= RC;
            for (size_t v=0; v<V; v++) { Y[v] =  3.0 + (Y[v]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
        }
        free(x1); free(xni);
    }
    else
    {
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        double sm1, sm2, sm4, *x1;
        if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t g=0; g<G; g++, X+=B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=J, Y++)
            {
                cblas_dcopy((int)L,X,(int)K,x1,1);
                sm1 = cblas_ddot((int)L,x1,1,&o,0);
                cblas_daxpy((int)L,-sm1,&ni,0,x1,1);
                sm2 = sm4 = 0.0;
                for (size_t l=0; l<L; l++) { double x = x1[l]*x1[l]; sm2 += x; sm4 += x*x; }
                *Y = L * sm4 / (sm2*sm2);
            }
        }
        if (!biased)
        {
            Y -= G*B;
            for (size_t n=0; n<N/L; n++) { Y[n] =  3.0 + (Y[n]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
        }
        free(x1);
    }

    return 0;
}


int kurtosis_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in kurtosis_c: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ni[2] = {-1.0f/L,0.0f};

    if (L<4)  { fprintf(stderr,"error in kurtosis_c: L must be > 3 \n"); return 1; }
    else if (L==N)
    {
        float sm1[2] = {0.0f,0.0f}, sm2[2] = {0.0f,0.0f}, sm4[2] = {0.0f,0.0f};
        float tmp, *x2;
        if (!(x2=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in kurtosis_c: problem with malloc. "); perror("malloc"); return 1; }
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
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        float tmp, sm1[2], sm2[2], sm4[2], *x1, *x2;
        if (!(x1=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in kurtosis_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(x2=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in kurtosis_c: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t g=0; g<G; g++, X+=2*B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=2*J)
            {
                cblas_ccopy((int)L,X,(int)K,x1,1);
                sm1[0] = sm1[1] = sm2[0] = sm2[1] = sm4[0] = sm4[1] = 0.0f;
                for (size_t l=0; l<2*L; l+=2) { sm1[0] += x1[l]; sm1[1] += x1[l+1]; }
                cblas_caxpy((int)L,sm1,ni,0,x1,1);
                cblas_cdotc_sub((int)L,x1,1,x1,1,(_Complex float *)sm2);
                for (size_t l=0; l<2*L; l+=2)
                {
                    x2[l] = x1[l]*x1[l] - x1[l+1]*x1[l+1];
                    x2[l+1] = x1[l]*x1[l+1] + x1[l+1]*x1[l];
                    tmp = x2[l]*x1[l] - x2[l+1]*x1[l+1];
                    x2[l+1] = x2[l]*x1[l+1] + x2[l+1]*x1[l];
                    x2[l] = tmp;
                }
                cblas_cdotu_sub((int)L,x2,1,x1,1,(_Complex float *)sm4);
                *Y++ = L*sm4[0]/(sm2[0]*sm2[0]);
                *Y++ = L*sm4[1]/(sm2[0]*sm2[0]);
            }
        }
        if (!biased)
        {
            Y -= 2*G*B;
            for (size_t v=0; v<2*N/L; v+=2)
            {
                Y[v] = 3.0f + (Y[v]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
                Y[v+1] *= (L+1)*(L-1) / (float)((L-2)*(L-3));
            }
        }
        free(x1); free(x2);
    }

    return 0;
}


int kurtosis_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in kurtosis_z: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni[2] = {-1.0/L,0.0};

    if (L<4)  { fprintf(stderr,"error in kurtosis_z: L must be > 3 \n"); return 1; }
    else if (L==N)
    {
        double sm1[2] = {0.0,0.0}, sm2[2] = {0.0,0.0}, sm4[2] = {0.0,0.0};
        double tmp, *x2;
        if (!(x2=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in kurtosis_z: problem with malloc. "); perror("malloc"); return 1; }
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
        const size_t B = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t G = N / (B*L);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        double tmp, sm1[2], sm2[2], sm4[2], *x1, *x2;
        if (!(x1=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in kurtosis_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(x2=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in kurtosis_z: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t g=0; g<G; g++, X+=2*B*(L-J))
        {
            for (size_t b=0; b<B; b++, X+=2*J)
            {
                cblas_zcopy((int)L,X,(int)K,x1,1);
                sm1[0] = sm1[1] = sm2[0] = sm2[1] = sm4[0] = sm4[1] = 0.0;
                for (size_t l=0; l<2*L; l+=2) { sm1[0] += x1[l]; sm1[1] += x1[l+1]; }
                cblas_zaxpy((int)L,sm1,ni,0,x1,1);
                cblas_zdotc_sub((int)L,x1,1,x1,1,(_Complex double *)sm2);
                for (size_t l=0; l<2*L; l+=2)
                {
                    x2[l] = x1[l]*x1[l] - x1[l+1]*x1[l+1];
                    x2[l+1] = x1[l]*x1[l+1] + x1[l+1]*x1[l];
                    tmp = x2[l]*x1[l] - x2[l+1]*x1[l+1];
                    x2[l+1] = x2[l]*x1[l+1] + x2[l+1]*x1[l];
                    x2[l] = tmp;
                }
                cblas_zdotu_sub((int)L,x2,1,x1,1,(_Complex double *)sm4);
                *Y++ = L*sm4[0]/(sm2[0]*sm2[0]);
                *Y++ = L*sm4[1]/(sm2[0]*sm2[0]);
            }
        }
        if (!biased)
        {
            Y -= 2*G*B;
            for (size_t v=0; v<2*N/L; v+=2)
            {
                Y[v] = 3.0 + (Y[v]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
                Y[v+1] *= (L+1)*(L-1) / (double)((L-2)*(L-3));
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
