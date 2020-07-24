//Vec2scalar (reduction) operation.
//Gets skewness for each vector in X along dim.
//This works in place.

//For complex case, output is complex.
//I follow the Octave convention for complex skewness (but see literature for other ideas later).

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int skewness_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int skewness_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int skewness_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int skewness_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);


int skewness_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in skewness_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float o = 1.0f, ni = 1.0f/L;
    const float w = (biased) ? sqrtf(L) : L*sqrtf(L-1)/(L-2);
    float x2, sm2 = 0.0f, sm3 = 0.0f;

    if (N==0) {}
    else if (L<3) { fprintf(stderr,"error in skewness_s: L must be > 2\n"); return 1; }
    else if (L==N)
    {
        *Y = cblas_sdot((int)L,X,1,&o,0) * ni;
        for (size_t l=0; l<L; l++) { *X -= *Y; x2 = *X**X; sm2 += x2; sm3 += x2**X++; }
        *Y = w * sm3 / (sm2*sqrtf(sm2));
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
                for (size_t v=0; v<V; v++)
                {
                    *Y = cblas_sdot((int)L,X,1,&o,0) * ni;
                    sm2 = sm3 = 0.0f;
                    for (size_t l=0; l<L; l++) { *X -= *Y; x2 = *X**X; sm2 += x2; sm3 += x2**X++; }
                    *Y++ = w * sm3 / (sm2*sqrtf(sm2));
                }
            }
            else
            {
                float *x1;
                if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with malloc. "); perror("malloc"); return 1; }
                cblas_scopy((int)L,&o,0,x1,1);
                cblas_sgemv(CblasColMajor,CblasTrans,(int)L,(int)V,o,X,(int)L,x1,1,0.0f,Y,1);
                for (size_t v=0; v<V; v++)
                {
                    *Y *= ni;
                    sm2 = sm3 = 0.0f;
                    for (size_t l=0; l<L; l++) { *X -= *Y; x2 = *X**X; sm2 += x2; sm3 += x2**X++; }
                    *Y++ = w * sm3 / (sm2*sqrtf(sm2));
                }
                free(x1);
            }
        }
        else
        {
            float *x1;
            if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, x1-=L, X++)
                {
                    cblas_scopy((int)L,X,(int)K,x1,1);
                    *Y = cblas_sdot((int)L,x1,1,&o,0) * ni;
                    sm2 = sm3 = 0.0f;
                    for (size_t l=0; l<L; l++) { *x1 -= *Y; x2 = *x1**x1; sm2 += x2; sm3 += x2**x1++; }
                    *Y++ = w * sm3 / (sm2*sqrtf(sm2));
                }
            }
            free(x1);
        }
    }

    return 0;
}


int skewness_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in skewness_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double o = 1.0, ni = 1.0/L;
    const double w = (biased) ? sqrt(L) : L*sqrt(L-1)/(L-2);
    double x2, sm2 = 0.0, sm3 = 0.0;

    if (N==0) {}
    else if (L<3) { fprintf(stderr,"error in skewness_d: L must be > 2\n"); return 1; }
    else if (L==N)
    {
        *Y = cblas_ddot((int)L,X,1,&o,0) * ni;
        for (size_t l=0; l<L; l++) { *X -= *Y; x2 = *X**X; sm2 += x2; sm3 += x2**X++; }
        *Y = w * sm3 / (sm2*sqrt(sm2));
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
                for (size_t v=0; v<V; v++)
                {
                    *Y = cblas_ddot((int)L,X,1,&o,0) * ni;
                    sm2 = sm3 = 0.0;
                    for (size_t l=0; l<L; l++) { *X -= *Y; x2 = *X**X; sm2 += x2; sm3 += x2**X++; }
                    *Y++ = w * sm3 / (sm2*sqrt(sm2));
                }
            }
            else
            {
                double *x1;
                if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with malloc. "); perror("malloc"); return 1; }
                cblas_dcopy((int)L,&o,0,x1,1);
                cblas_dgemv(CblasColMajor,CblasTrans,(int)L,(int)V,o,X,(int)L,x1,1,0.0,Y,1);
                for (size_t v=0; v<V; v++)
                {
                    *Y *= ni;
                    sm2 = sm3 = 0.0;
                    for (size_t l=0; l<L; l++) { *X -= *Y; x2 = *X**X; sm2 += x2; sm3 += x2**X++; }
                    *Y++ = w * sm3 / (sm2*sqrt(sm2));
                }
                free(x1);
            }
        }
        else
        {
            double *x1;
            if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, x1-=L, X++)
                {
                    cblas_dcopy((int)L,X,(int)K,x1,1);
                    *Y = cblas_ddot((int)L,x1,1,&o,0) * ni;
                    sm2 = sm3 = 0.0;
                    for (size_t l=0; l<L; l++) { *x1 -= *Y; x2 = *x1**x1; sm2 += x2; sm3 += x2**x1++; }
                    *Y++ = w * sm3 / (sm2*sqrt(sm2));
                }
            }
            free(x1);
        }
    }

    return 0;
}


int skewness_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in skewness_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ni[2] = {-1.0f/L,0.0f};
    const float w = (biased) ? sqrtf(L) : L*sqrtf(L-1)/(L-2);
    float x2[2], sm1[2] = {0.0f,0.0f}, sm2[2] = {0.0f,0.0f}, sm3[2] = {0.0f,0.0f};

    if (N==0) {}
    else if (L<3) { fprintf(stderr,"error in skewness_c: L must be > 2\n"); return 1; }
    else if (L==N)
    {
        for (size_t l=0; l<2*L; l+=2) { sm1[0] += X[l]; sm1[1] += X[l+1]; }
        cblas_caxpy((int)L,sm1,ni,0,X,1);
        cblas_cdotc_sub((int)L,X,1,X,1,(_Complex float *)sm2);
        for (size_t l=0; l<2*L; l+=2)
        {
            x2[0] = X[l]*X[l] - X[l+1]*X[l+1];
            x2[1] = X[l]*X[l+1] + X[l+1]*X[l];
            sm3[0] += x2[0]*X[l] - x2[1]*X[l+1];
            sm3[1] += x2[0]*X[l+1] + x2[1]*X[l];
        }
        Y[0] = w * sm3[0] / (sm2[0]*sqrtf(sm2[0]));
        Y[1] = w * sm3[1] / (sm2[0]*sqrtf(sm2[0]));
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++)
            {
                sm1[0] = sm1[1] = sm2[0] = sm2[1] = sm3[0] = sm3[1] = 0.0f;
                for (size_t l=0; l<2*L; l+=2) { sm1[0] += X[l]; sm1[1] += X[l+1]; }
                cblas_caxpy((int)L,sm1,ni,0,X,1);
                cblas_cdotc_sub((int)L,X,1,X,1,(_Complex float *)sm2);
                for (size_t l=0; l<2*L; l+=2)
                {
                    x2[0] = X[l]*X[l] - X[l+1]*X[l+1];
                    x2[1] = X[l]*X[l+1] + X[l+1]*X[l];
                    sm3[0] += x2[0]*X[l] - x2[1]*X[l+1];
                    sm3[1] += x2[0]*X[l+1] + x2[1]*X[l];
                }
                *Y++ = w * sm3[0] / (sm2[0]*sqrtf(sm2[0]));
                *Y++ = w * sm3[1] / (sm2[0]*sqrtf(sm2[0]));
            }
        }
        else
        {
            float *x1;
            if (!(x1=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in skewness_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X+=2)
                {
                    cblas_ccopy((int)L,X,(int)K,x1,1);
                    sm1[0] = sm1[1] = sm2[0] = sm2[1] = sm3[0] = sm3[1] = 0.0f;
                    for (size_t l=0; l<2*L; l+=2) { sm1[0] += x1[l]; sm1[1] += x1[l+1]; }
                    cblas_caxpy((int)L,sm1,ni,0,x1,1);
                    cblas_cdotc_sub((int)L,x1,1,x1,1,(_Complex float *)sm2);
                    for (size_t l=0; l<2*L; l+=2)
                    {
                        x2[0] = x1[l]*x1[l] - x1[l+1]*x1[l+1];
                        x2[1] = x1[l]*x1[l+1] + x1[l+1]*x1[l];
                        sm3[0] += x2[0]*x1[l] - x2[1]*x1[l+1];
                        sm3[1] += x2[0]*x1[l+1] + x2[1]*x1[l];
                    }
                    *Y++ = w * sm3[0] / (sm2[0]*sqrtf(sm2[0]));
                    *Y++ = w * sm3[1] / (sm2[0]*sqrtf(sm2[0]));
                }
            }
            free(x1);
        }
    }

    return 0;
}


int skewness_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in skewness_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni[2] = {-1.0/L,0.0};
    const double w = (biased) ? sqrt(L) : L*sqrt(L-1)/(L-2);
    double x2[2], sm1[2] = {0.0,0.0}, sm2[2] = {0.0,0.0}, sm3[2] = {0.0,0.0};

    if (N==0) {}
    else if (L<3) { fprintf(stderr,"error in skewness_z: L must be > 2\n"); return 1; }
    else if (L==N)
    {
        for (size_t l=0; l<2*L; l+=2) { sm1[0] += X[l]; sm1[1] += X[l+1]; }
        cblas_zaxpy((int)L,sm1,ni,0,X,1);
        cblas_zdotc_sub((int)L,X,1,X,1,(_Complex double *)sm2);
        for (size_t l=0; l<2*L; l+=2)
        {
            x2[0] = X[l]*X[l] - X[l+1]*X[l+1];
            x2[1] = X[l]*X[l+1] + X[l+1]*X[l];
            sm3[0] += x2[0]*X[l] - x2[1]*X[l+1];
            sm3[1] += x2[0]*X[l+1] + x2[1]*X[l];
        }
        Y[0] = w * sm3[0] / (sm2[0]*sqrt(sm2[0]));
        Y[1] = w * sm3[1] / (sm2[0]*sqrt(sm2[0]));
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; v++)
            {
                sm1[0] = sm1[1] = sm2[0] = sm2[1] = sm3[0] = sm3[1] = 0.0;
                for (size_t l=0; l<2*L; l+=2) { sm1[0] += X[l]; sm1[1] += X[l+1]; }
                cblas_zaxpy((int)L,sm1,ni,0,X,1);
                cblas_zdotc_sub((int)L,X,1,X,1,(_Complex double *)sm2);
                for (size_t l=0; l<2*L; l+=2)
                {
                    x2[0] = X[l]*X[l] - X[l+1]*X[l+1];
                    x2[1] = X[l]*X[l+1] + X[l+1]*X[l];
                    sm3[0] += x2[0]*X[l] - x2[1]*X[l+1];
                    sm3[1] += x2[0]*X[l+1] + x2[1]*X[l];
                }
                *Y++ = w * sm3[0] / (sm2[0]*sqrt(sm2[0]));
                *Y++ = w * sm3[1] / (sm2[0]*sqrt(sm2[0]));
            }
        }
        else
        {
            double *x1;
            if (!(x1=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in skewness_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X+=2)
                {
                    cblas_zcopy((int)L,X,(int)K,x1,1);
                    sm1[0] = sm1[1] = sm2[0] = sm2[1] = sm3[0] = sm3[1] = 0.0;
                    for (size_t l=0; l<2*L; l+=2) { sm1[0] += x1[l]; sm1[1] += x1[l+1]; }
                    cblas_zaxpy((int)L,sm1,ni,0,x1,1);
                    cblas_zdotc_sub((int)L,x1,1,x1,1,(_Complex double *)sm2);
                    for (size_t l=0; l<2*L; l+=2)
                    {
                        x2[0] = x1[l]*x1[l] - x1[l+1]*x1[l+1];
                        x2[1] = x1[l]*x1[l+1] + x1[l+1]*x1[l];
                        sm3[0] += x2[0]*x1[l] - x2[1]*x1[l+1];
                        sm3[1] += x2[0]*x1[l+1] + x2[1]*x1[l];
                    }
                    *Y++ = w * sm3[0] / (sm2[0]*sqrt(sm2[0]));
                    *Y++ = w * sm3[1] / (sm2[0]*sqrt(sm2[0]));
                }
            }
            free(x1);
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
