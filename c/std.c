//Vec2scalar (reduction) operation.
//Gets standard deviation of each row or col of X according to dim.
//This works in place.

//For complex case, output is real.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int std_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int std_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int std_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int std_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);


int std_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in std_s: dim must be in [0 3]\n"); return 1; }

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
        *Y = cblas_snrm2((int)L,X,1) * den;
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
                for (size_t v=0; v<V; v++, X+=L)
                {
                    *Y = cblas_sdot((int)L,X,1,&o,0) * ni;
                    cblas_saxpy((int)L,-*Y,&o,0,X,1);
                    *Y++ = cblas_snrm2((int)L,X,1) * den;
                }
            }
            else
            {
                float *x1, *xni;
                if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in std_s: problem with malloc. "); perror("malloc"); return 1; }
                if (!(xni=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in std_s: problem with malloc. "); perror("malloc"); return 1; }
                cblas_scopy((int)L,&o,0,x1,1); cblas_scopy((int)L,&ni,0,xni,1);
                cblas_sgemv(CblasColMajor,CblasTrans,(int)L,(int)V,o,X,(int)L,x1,1,0.0f,Y,1);
                cblas_sgemm(CblasColMajor,CblasTrans,CblasTrans,(int)L,(int)V,1,-o,xni,1,Y,(int)V,o,X,(int)L);
                for (size_t v=0; v<V; v++, X+=L) { *Y++ = cblas_snrm2((int)L,X,1) * den; }
                free(x1); free(xni);
            }
        }
        else if (G==1)
        {
            float *x1, *xni;
            if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in std_s: problem with malloc. "); perror("malloc"); return 1; }
            if (!(xni=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in std_s: problem with malloc. "); perror("malloc"); return 1; }
            cblas_scopy((int)L,&o,0,x1,1); cblas_scopy((int)L,&ni,0,xni,1);
            cblas_sgemv(CblasRowMajor,CblasTrans,(int)L,(int)V,o,X,(int)V,x1,1,z,Y,1);
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)L,(int)V,1,-o,xni,1,Y,(int)V,o,X,(int)V);
            for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
            cblas_sgemv(CblasRowMajor,CblasTrans,(int)L,(int)V,o,X,(int)V,x1,1,z,Y,1);
            for (size_t v=0; v<V; v++) { Y[v] = sqrtf(Y[v]) * den; }
            free(x1); free(xni);
        }
        else
        {
            float *x1;
            if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in std_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, X++)
                {
                    cblas_scopy((int)L,X,(int)K,x1,1);
                    *Y = cblas_sdot((int)L,x1,1,&o,0) * ni;
                    cblas_saxpy((int)L,-*Y,&o,0,x1,1);
                    *Y++ = cblas_snrm2((int)L,x1,1) * den;
                }
            }
            free(x1);
        }
    }
    
    return 0;
}


int std_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in std_d: dim must be in [0 3]\n"); return 1; }

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
        *Y = cblas_dnrm2((int)L,X,1) * den;
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
                for (size_t v=0; v<V; v++, X+=L)
                {
                    *Y = cblas_ddot((int)L,X,1,&o,0) * ni;
                    cblas_daxpy((int)L,-*Y,&o,0,X,1);
                    *Y++ = cblas_dnrm2((int)L,X,1) * den;
                }
            }
            else
            {
                double *x1, *xni;
                if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in std_d: problem with malloc. "); perror("malloc"); return 1; }
                if (!(xni=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in std_d: problem with malloc. "); perror("malloc"); return 1; }
                cblas_dcopy((int)L,&o,0,x1,1); cblas_dcopy((int)L,&ni,0,xni,1);
                cblas_dgemv(CblasColMajor,CblasTrans,(int)L,(int)V,o,X,(int)L,x1,1,z,Y,1);
                cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,(int)L,(int)V,1,-o,xni,1,Y,(int)V,o,X,(int)L);
                for (size_t v=0; v<V; v++, X+=L) { *Y++ = cblas_dnrm2((int)L,X,1) * den; }
                free(x1); free(xni);
            }
        }
        else if (G==1)
        {
            double *x1, *xni;
            if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in std_d: problem with malloc. "); perror("malloc"); return 1; }
            if (!(xni=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in std_d: problem with malloc. "); perror("malloc"); return 1; }
            cblas_dcopy((int)L,&o,0,x1,1); cblas_dcopy((int)L,&ni,0,xni,1);
            cblas_dgemv(CblasRowMajor,CblasTrans,(int)L,(int)V,o,X,(int)V,x1,1,z,Y,1);
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)L,(int)V,1,-o,xni,1,Y,(int)V,o,X,(int)V);
            for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
            cblas_dgemv(CblasRowMajor,CblasTrans,(int)L,(int)V,o,X,(int)V,x1,1,z,Y,1);
            for (size_t v=0; v<V; v++) { Y[v] = sqrt(Y[v]) * den; }
            free(x1); free(xni);
        }
        else
        {
            double *x1;
            if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in std_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, X++)
                {
                    cblas_dcopy((int)L,X,(int)K,x1,1);
                    *Y = cblas_ddot((int)L,x1,1,&o,0) * ni;
                    cblas_daxpy((int)L,-*Y,&o,0,x1,1);
                    *Y++ = cblas_dnrm2((int)L,x1,1) * den;
                }
            }
            free(x1);
        }
    }
    
    return 0;
}


int std_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in std_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = (biased) ? 1.0f/sqrtf(L) : 1.0f/sqrtf(L-1);
    const float ni[2] = {-1.0f/L,0.0f};

    if (N==0) {}
    else if (L==1)
    {
        const float z[2] = {0.0f,0.0f};
        cblas_ccopy((int)N,z,0,Y,1);
    }
    else if (L==N)
    {
        float mn[2] = {0.0f,0.0f};
        for (size_t n=0; n<2*L; n+=2) { mn[0] += X[n]; mn[1] += X[n+1]; }
        cblas_caxpy((int)L,mn,ni,0,X,1);
        *Y = cblas_scnrm2((int)L,X,1) * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            float mn[2], *xni;
            if (!(xni=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in std_c: problem with malloc. "); perror("malloc"); return 1; }
            cblas_ccopy((int)L,ni,0,xni,1);
            for (size_t v=0; v<V; v++, X+=2*L)
            {
                mn[0] = mn[1] = 0.0f;
                for (size_t l=0; l<L; l++) { mn[0] += *X++; mn[1] += *X++; }
                X -= 2*L;
                cblas_caxpy((int)L,mn,ni,0,X,1);
                *Y++ = cblas_scnrm2((int)L,X,1) * den;
            }
            free(xni);
        }
        else
        {
            float mn[2], *x1;
            if (!(x1=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in std_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X+=2)
                {
                    cblas_ccopy((int)L,X,(int)K,x1,1);
                    mn[0] = mn[1] = 0.0f;
                    for (size_t l=0; l<2*L; l+=2) { mn[0] += x1[l]; mn[1] += x1[l+1]; }
                    cblas_caxpy((int)L,mn,ni,0,x1,1);
                    *Y++ = cblas_scnrm2((int)L,x1,1) * den;
                }
            }
            free(x1);
        }
    }
    
    return 0;
}


int std_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in std_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = (biased) ? 1.0/sqrt(L) : 1.0/sqrt(L-1);
    const double ni[2] = {-1.0/L,0.0};

    if (N==0) {}
    else if (L==1)
    {
        const double z[2] = {0.0,0.0};
        cblas_zcopy((int)N,z,0,Y,1);
    }
    else if (L==N)
    {
        double mn[2] = {0.0,0.0};
        for (size_t n=0; n<2*L; n+=2) { mn[0] += X[n]; mn[1] += X[n+1]; }
        cblas_zaxpy((int)L,mn,ni,0,X,1);
        *Y = cblas_dznrm2((int)L,X,1) * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            double mn[2], *xni;
            if (!(xni=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in std_z: problem with malloc. "); perror("malloc"); return 1; }
            cblas_zcopy((int)L,ni,0,xni,1);
            for (size_t v=0; v<V; v++, X+=2*L)
            {
                mn[0] = mn[1] = 0.0;
                for (size_t l=0; l<L; l++) { mn[0] += *X++; mn[1] += *X++; }
                X -= 2*L;
                cblas_zaxpy((int)L,mn,ni,0,X,1);
                *Y++ = cblas_dznrm2((int)L,X,1);
            }
            free(xni);
        }
        else
        {
            double mn[2], *x1;
            if (!(x1=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in std_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X+=2)
                {
                    cblas_zcopy((int)L,X,(int)K,x1,1);
                    mn[0] = mn[1] = 0.0;
                    for (size_t l=0; l<2*L; l+=2) { mn[0] += x1[l]; mn[1] += x1[l+1]; }
                    cblas_zaxpy((int)L,mn,ni,0,x1,1);
                    *Y++ = cblas_dznrm2((int)L,x1,1) * den;
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
