//Vec2scalar (reduction) operation.
//Gets the kurtosis for each vector in X along dim.
//This works in place.

//For complex case, output is complex.
//I follow the Octave convention for complex kurtosis (but see literature for other ideas later).

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

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float o = 1.0f, ni = 1.0f/L;
    float x2, sm2 = 0.0f, sm4 = 0.0f;

    if (N==0) {}
    else if (L<4) { fprintf(stderr,"error in kurtosis_s: L must be > 3\n"); return 1; }
    else if (L==N)
    {
        *Y = cblas_sdot((int)L,X,1,&o,0) * ni;
        for (size_t l=0; l<L; l++, X++) { *X -= *Y; x2 = *X**X; sm2 += x2; sm4 += x2*x2; }
        *Y = L * sm4 / (sm2*sm2);
        if (!biased) { *Y =  3.0f + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            if (V<250)
            {
                for (size_t v=0; v<V; v++, Y++)
                {
                    *Y = cblas_sdot((int)L,X,1,&o,0) * ni;
                    sm2 = sm4 = 0.0f;
                    for (size_t l=0; l<L; l++, X++) { *X -= *Y; x2 = *X**X; sm2 += x2; sm4 += x2*x2; }
                    *Y = L * sm4 / (sm2*sm2);
                    if (!biased) { *Y =  3.0f + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
                }
            }
            else
            {
                float *x1;
                if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with malloc. "); perror("malloc"); return 1; }
                cblas_scopy((int)L,&o,0,x1,1);
                cblas_sgemv(CblasColMajor,CblasTrans,(int)L,(int)V,o,X,(int)L,x1,1,0.0f,Y,1);
                for (size_t v=0; v<V; v++, Y++)
                {
                    *Y *= ni;
                    sm2 = sm4 = 0.0f;
                    for (size_t l=0; l<L; l++, X++) { *X -= *Y; x2 = *X**X; sm2 += x2; sm4 += x2*x2; }
                    *Y = L * sm4 / (sm2*sm2);
                    if (!biased) { *Y =  3.0f + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
                }
                free(x1);
            }
        }
        else
        {
            float *x1;
            if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, x1-=L, X++, Y++)
                {
                    cblas_scopy((int)L,X,(int)K,x1,1);
                    *Y = cblas_sdot((int)L,x1,1,&o,0) * ni;
                    sm2 = sm4 = 0.0f;
                    for (size_t l=0; l<L; l++, x1++) { *x1 -= *Y; x2 = *x1**x1; sm2 += x2; sm4 += x2*x2; }
                    *Y = L * sm4 / (sm2*sm2);
                    if (!biased) { *Y =  3.0f + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
                }
            }
            free(x1);
        }
    }

    return 0;
}


int kurtosis_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in kurtosis_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double o = 1.0, ni = 1.0/L;
    double x2, sm2 = 0.0, sm4 = 0.0;

    if (N==0) {}
    else if (L<4) { fprintf(stderr,"error in kurtosis_d: L must be > 3\n"); return 1; }
    else if (L==N)
    {
        *Y = cblas_ddot((int)L,X,1,&o,0) * ni;
        for (size_t l=0; l<L; l++, X++) { *X -= *Y; x2 = *X**X; sm2 += x2; sm4 += x2*x2; }
        *Y = L * sm4 / (sm2*sm2);
        if (!biased) { *Y =  3.0 + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            if (V<250)
            {
                for (size_t v=0; v<V; v++, Y++)
                {
                    *Y = cblas_ddot((int)L,X,1,&o,0) * ni;
                    sm2 = sm4 = 0.0;
                    for (size_t l=0; l<L; l++, X++) { *X -= *Y; x2 = *X**X; sm2 += x2; sm4 += x2*x2; }
                    *Y = L * sm4 / (sm2*sm2);
                    if (!biased) { *Y =  3.0 + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
                }
            }
            else
            {
                double *x1;
                if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with malloc. "); perror("malloc"); return 1; }
                cblas_dcopy((int)L,&o,0,x1,1);
                cblas_dgemv(CblasColMajor,CblasTrans,(int)L,(int)V,o,X,(int)L,x1,1,0.0,Y,1);
                for (size_t v=0; v<V; v++, Y++)
                {
                    *Y *= ni;
                    sm2 = sm4 = 0.0;
                    for (size_t l=0; l<L; l++, X++) { *X -= *Y; x2 = *X**X; sm2 += x2; sm4 += x2*x2; }
                    *Y = L * sm4 / (sm2*sm2);
                    if (!biased) { *Y =  3.0 + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
                }
                free(x1);
            }
        }
        else
        {
            double *x1;
            if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, x1-=L, X++, Y++)
                {
                    cblas_dcopy((int)L,X,(int)K,x1,1);
                    *Y = cblas_ddot((int)L,x1,1,&o,0) * ni;
                    sm2 = sm4 = 0.0;
                    for (size_t l=0; l<L; l++, x1++) { *x1 -= *Y; x2 = *x1**x1; sm2 += x2; sm4 += x2*x2; }
                    *Y = L * sm4 / (sm2*sm2);
                    if (!biased) { *Y =  3.0 + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
                }
            }
            free(x1);
        }
    }

    return 0;
}


int kurtosis_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in kurtosis_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ni[2] = {-1.0f/L,0.0f};
    float sm1[2] = {0.0f,0.0f}, sm2[2] = {0.0f,0.0f}, sm4[2] = {0.0f,0.0f};
    float tmp, *x2;
    if (!(x2=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in kurtosis_c: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L<4) { fprintf(stderr,"error in kurtosis_c: L must be > 3\n"); return 1; }
    else if (L==N)
    {
        for (size_t l=0; l<2*L; l+=2) { sm1[0] += X[l]; sm1[1] += X[l+1]; }
        cblas_caxpy((int)L,sm1,ni,0,X,1);
        cblas_cdotc_sub((int)L,X,1,X,1,(_Complex float *)sm2);
        for (size_t l=0; l<2*L; l+=2)
        {
            x2[l] = X[l]*X[l] - X[l+1]*X[l+1];
            x2[l+1] = X[l]*X[l+1] + X[l+1]*X[l];
            tmp = x2[l]*X[l] - x2[l+1]*X[l+1];
            x2[l+1] = x2[l]*X[l+1] + x2[l+1]*X[l];
            x2[l] = tmp;
        }
        cblas_cdotu_sub((int)L,x2,1,X,1,(_Complex float *)sm4);
        Y[0] = L * sm4[0] / (sm2[0]*sm2[0]);
        Y[1] = L * sm4[1] / (sm2[0]*sm2[0]);
        if (!biased)
        {
            Y[0] = 3.0f + (Y[0]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
            Y[1] *= (L+1)*(L-1) / (float)((L-2)*(L-3));
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, X+=2*L, Y+=2)
            {
                sm1[0] = sm1[1] = sm2[0] = sm2[1] = sm4[0] = sm4[1] = 0.0f;
                for (size_t l=0; l<2*L; l+=2) { sm1[0] += X[l]; sm1[1] += X[l+1]; }
                cblas_caxpy((int)L,sm1,ni,0,X,1);
                cblas_cdotc_sub((int)L,X,1,X,1,(_Complex float *)sm2);
                for (size_t l=0; l<2*L; l+=2)
                {
                    x2[l] = X[l]*X[l] - X[l+1]*X[l+1];
                    x2[l+1] = X[l]*X[l+1] + X[l+1]*X[l];
                    tmp = x2[l]*X[l] - x2[l+1]*X[l+1];
                    x2[l+1] = x2[l]*X[l+1] + x2[l+1]*X[l];
                    x2[l] = tmp;
                }
                cblas_cdotu_sub((int)L,x2,1,X,1,(_Complex float *)sm4);
                Y[0] = L * sm4[0] / (sm2[0]*sm2[0]);
                Y[1] = L * sm4[1] / (sm2[0]*sm2[0]);
                if (!biased)
                {
                    Y[0] = 3.0f + (Y[0]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
                    Y[1] *= (L+1)*(L-1) / (float)((L-2)*(L-3));
                }
            }
        }
        else
        {
            float *x1;
            if (!(x1=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in kurtosis_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X+=2, Y+=2)
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
                    Y[0] = L * sm4[0] / (sm2[0]*sm2[0]);
                    Y[1] = L * sm4[1] / (sm2[0]*sm2[0]);
                    if (!biased)
                    {
                        Y[0] = 3.0f + (Y[0]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
                        Y[1] *= (L+1)*(L-1) / (float)((L-2)*(L-3));
                    }
                }
            }
            free(x1);
        }
    }

    free(x2);
    return 0;
}


int kurtosis_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in kurtosis_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni[2] = {-1.0/L,0.0};
    double sm1[2] = {0.0,0.0}, sm2[2] = {0.0,0.0}, sm4[2] = {0.0,0.0};
    double tmp, *x2;
    if (!(x2=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in kurtosis_z: problem with malloc. "); perror("malloc"); return 1; }

    if (N==0) {}
    else if (L<4) { fprintf(stderr,"error in kurtosis_z: L must be > 3\n"); return 1; }
    else if (L==N)
    {
        for (size_t l=0; l<2*L; l+=2) { sm1[0] += X[l]; sm1[1] += X[l+1]; }
        cblas_zaxpy((int)L,sm1,ni,0,X,1);
        cblas_zdotc_sub((int)L,X,1,X,1,(_Complex double *)sm2);
        for (size_t l=0; l<2*L; l+=2)
        {
            x2[l] = X[l]*X[l] - X[l+1]*X[l+1];
            x2[l+1] = X[l]*X[l+1] + X[l+1]*X[l];
            tmp = x2[l]*X[l] - x2[l+1]*X[l+1];
            x2[l+1] = x2[l]*X[l+1] + x2[l+1]*X[l];
            x2[l] = tmp;
        }
        cblas_zdotu_sub((int)L,x2,1,X,1,(_Complex double *)sm4);
        Y[0] = L * sm4[0] / (sm2[0]*sm2[0]);
        Y[1] = L * sm4[1] / (sm2[0]*sm2[0]);
        if (!biased)
        {
            Y[0] = 3.0 + (Y[0]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
            Y[1] *= (L+1)*(L-1) / (double)((L-2)*(L-3));
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++, X+=2*L, Y+=2)
            {
                sm1[0] = sm1[1] = sm2[0] = sm2[1] = sm4[0] = sm4[1] = 0.0;
                for (size_t l=0; l<2*L; l+=2) { sm1[0] += X[l]; sm1[1] += X[l+1]; }
                cblas_zaxpy((int)L,sm1,ni,0,X,1);
                cblas_zdotc_sub((int)L,X,1,X,1,(_Complex double *)sm2);
                for (size_t l=0; l<2*L; l+=2)
                {
                    x2[l] = X[l]*X[l] - X[l+1]*X[l+1];
                    x2[l+1] = X[l]*X[l+1] + X[l+1]*X[l];
                    tmp = x2[l]*X[l] - x2[l+1]*X[l+1];
                    x2[l+1] = x2[l]*X[l+1] + x2[l+1]*X[l];
                    x2[l] = tmp;
                }
                cblas_zdotu_sub((int)L,x2,1,X,1,(_Complex double *)sm4);
                Y[0] = L * sm4[0] / (sm2[0]*sm2[0]);
                Y[1] = L * sm4[1] / (sm2[0]*sm2[0]);
                if (!biased)
                {
                    Y[0] = 3.0 + (Y[0]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
                    Y[1] *= (L+1)*(L-1) / (double)((L-2)*(L-3));
                }
            }
        }
        else
        {
            double *x1;
            if (!(x1=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in kurtosis_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X+=2, Y+=2)
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
                    Y[0] = L * sm4[0] / (sm2[0]*sm2[0]);
                    Y[1] = L * sm4[1] / (sm2[0]*sm2[0]);
                    if (!biased)
                    {
                        Y[0] = 3.0 + (Y[0]*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
                        Y[1] *= (L+1)*(L-1) / (double)((L-2)*(L-3));
                    }
                }
            }
            free(x1);
        }
    }

    free(x2);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
