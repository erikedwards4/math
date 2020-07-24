//Vec2scalar (reduction) operation.
//Gets the coefficient of variation (std/mean) for each vector in X along dim.
//This works in place.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int coeff_var_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int coeff_var_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);


int coeff_var_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in coeff_var_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float z = 0.0f, o = 1.0f, ni = 1.0f/L;
    const float den = (biased) ? 1.0f/sqrtf(L) : 1.0f/sqrtf(L-1);

    if (N==0) {}
    else if (L==1) { cblas_scopy((int)N,&z,0,Y,1); }
    else if (L==N)
    {
        *Y = cblas_sdot((int)L,X,1,&o,0) * ni;
        cblas_saxpy((int)L,-*Y,&o,0,X,1);
        *Y = cblas_snrm2((int)L,X,1) * den / *Y;
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
                for (size_t v=0; v<V; v++, X+=L, Y++)
                {
                    *Y = cblas_sdot((int)L,X,1,&o,0) * ni;
                    cblas_saxpy((int)L,-*Y,&o,0,X,1);
                    *Y = cblas_snrm2((int)L,X,1) * den / *Y;
                }
            }
            else
            {
                float *x1;
                if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
                cblas_scopy((int)L,&o,0,x1,1);
                cblas_sgemv(CblasColMajor,CblasTrans,(int)L,(int)V,o,X,(int)L,x1,1,0.0f,Y,1);
                cblas_sscal((int)V,ni,Y,1);
                cblas_sgemm(CblasColMajor,CblasTrans,CblasTrans,(int)L,(int)V,1,-o,x1,1,Y,(int)V,o,X,(int)L);
                for (size_t v=0; v<V; v++, X+=L, Y++) { *Y = cblas_snrm2((int)L,X,1) * den / *Y; }
                free(x1);
            }
        }
        else if (G==1)
        {
            float *x1;
            if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
            cblas_scopy((int)L,&o,0,x1,1);
            cblas_sgemv(CblasRowMajor,CblasTrans,(int)L,(int)V,o,X,(int)V,x1,1,z,Y,1);
            cblas_sscal((int)V,ni,Y,1);
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)L,(int)V,1,-o,x1,1,Y,(int)V,o,X,(int)V);
            for (size_t v=0; v<V; v++, X++, Y++) { *Y = cblas_snrm2((int)L,X,(int)V) * den / *Y; }
            free(x1);
        }
        else
        {
            float *x1;
            if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, X++, Y++)
                {
                    cblas_scopy((int)L,X,(int)K,x1,1);
                    *Y = cblas_sdot((int)L,x1,1,&o,0) * ni;
                    cblas_saxpy((int)L,-*Y,&o,0,x1,1);
                    *Y = cblas_snrm2((int)L,x1,1) * den / *Y;
                }
            }
            free(x1);
        }
    }
    
    return 0;
}


int coeff_var_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in coeff_var_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double z = 0.0, o = 1.0, ni = 1.0/L;
    const double den = (biased) ? 1.0/sqrt(L) : 1.0/sqrt(L-1);

    if (N==0) {}
    else if (L==1) { cblas_dcopy((int)N,&z,0,Y,1); }
    else if (L==N)
    {
        *Y = cblas_ddot((int)L,X,1,&o,0) * ni;
        cblas_daxpy((int)L,-*Y,&o,0,X,1);
        *Y = cblas_dnrm2((int)L,X,1) * den / *Y;
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
                for (size_t v=0; v<V; v++, X+=L, Y++)
                {
                    *Y = cblas_ddot((int)L,X,1,&o,0) * ni;
                    cblas_daxpy((int)L,-*Y,&o,0,X,1);
                    *Y = cblas_dnrm2((int)L,X,1) * den / *Y;
                }
            }
            else
            {
                double *x1;
                if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
                cblas_dcopy((int)L,&o,0,x1,1);
                cblas_dgemv(CblasColMajor,CblasTrans,(int)L,(int)V,o,X,(int)L,x1,1,0.0,Y,1);
                cblas_dscal((int)V,ni,Y,1);
                cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,(int)L,(int)V,1,-o,x1,1,Y,(int)V,o,X,(int)L);
                for (size_t v=0; v<V; v++, X+=L, Y++) { *Y = cblas_dnrm2((int)L,X,1) * den / *Y; }
                free(x1);
            }
        }
        else if (G==1)
        {
            double *x1;
            if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
            cblas_dcopy((int)L,&o,0,x1,1);
            cblas_dgemv(CblasRowMajor,CblasTrans,(int)L,(int)V,o,X,(int)V,x1,1,z,Y,1);
            cblas_dscal((int)V,ni,Y,1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)L,(int)V,1,-o,x1,1,Y,(int)V,o,X,(int)V);
            for (size_t v=0; v<V; v++, X++, Y++) { *Y = cblas_dnrm2((int)L,X,(int)V) * den / *Y; }
            free(x1);
        }
        else
        {
            double *x1;
            if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, X++, Y++)
                {
                    cblas_dcopy((int)L,X,(int)K,x1,1);
                    *Y = cblas_ddot((int)L,x1,1,&o,0) * ni;
                    cblas_daxpy((int)L,-*Y,&o,0,x1,1);
                    *Y = cblas_dnrm2((int)L,x1,1) * den / *Y;
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
