//Vec2scalar (reduction) operation.
//Gets mean for each vector in X along dim.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int mean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float o = 1.0f, ni = 1.0f / L;

    if (N==0) {}
    else if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y = cblas_sdot((int)L,X,1,&o,0) * ni;
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
                for (size_t v=0; v<V; v++, X+=L) { *Y++ = cblas_sdot((int)L,X,1,&o,0) * ni; }
            }
            else
            {
                float *x1;
                if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in mean_s: problem with malloc. "); perror("malloc"); return 1; }
                cblas_scopy((int)L,&o,0,x1,1);
                cblas_sgemv(CblasColMajor,CblasTrans,(int)L,(int)V,o,X,(int)L,x1,1,0.0f,Y,1);
                cblas_sscal((int)V,ni,Y,1);
                free(x1);
            }
        }
        else if (G==1)
        {
            if (V<40)
            {
                for (size_t v=0; v<V; v++, X++) { *Y++ = cblas_sdot((int)L,X,(int)V,&o,0) * ni; }
            }
            else
            {
                float *x1;
                if (!(x1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in mean_s: problem with malloc. "); perror("malloc"); return 1; }
                cblas_scopy((int)L,&o,0,x1,1);
                cblas_sgemv(CblasRowMajor,CblasTrans,(int)L,(int)V,o,X,(int)V,x1,1,0.0f,Y,1);
                cblas_sscal((int)V,ni,Y,1);
                free(x1);
            }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, X++)
                {
                    *Y++ = cblas_sdot((int)L,X,(int)K,&o,0) * ni;
                }
            }
        }
    }

    return 0;
}


int mean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double o = 1.0, ni = 1.0 / L;

    if (N==0) {}
    else if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y = cblas_ddot((int)L,X,1,&o,0) * ni;
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
                for (size_t v=0; v<V; v++, X+=L) { *Y++ = cblas_ddot((int)L,X,1,&o,0) * ni; }
            }
            else
            {
                double *x1;
                if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in mean_d: problem with malloc. "); perror("malloc"); return 1; }
                cblas_dcopy((int)L,&o,0,x1,1);
                cblas_dgemv(CblasColMajor,CblasTrans,(int)L,(int)V,o,X,(int)L,x1,1,0.0,Y,1);
                cblas_dscal((int)V,ni,Y,1);
                free(x1);
            }
        }
        else if (G==1)
        {
            if (V<40)
            {
                for (size_t v=0; v<V; v++, X++) { *Y++ = cblas_ddot((int)L,X,(int)V,&o,0) * ni; }
            }
            else
            {
                double *x1;
                if (!(x1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in mean_d: problem with malloc. "); perror("malloc"); return 1; }
                cblas_dcopy((int)L,&o,0,x1,1);
                cblas_dgemv(CblasRowMajor,CblasTrans,(int)L,(int)V,o,X,(int)V,x1,1,0.0,Y,1);
                cblas_dscal((int)V,ni,Y,1);
                free(x1);
            }
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=B*(L-1))
            {
                for (size_t b=0; b<B; b++, X++)
                {
                    *Y++ = cblas_ddot((int)L,X,(int)K,&o,0) * ni;
                }
            }
        }
    }

    return 0;
}


int mean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ni = 1.0f / L;

    if (N==0) {}
    else if (L==1) { cblas_ccopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y++ = 0.0f; *Y-- = 0.0f;
        for (size_t l=0; l<L; l++) { *Y++ += *X++; *Y-- += *X++; }
        *Y++ *= ni; *Y-- *= ni;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};

        if (K==1 && (G==1 || B==1))
        {
            float *x1;
            if (!(x1=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in mean_c: problem with malloc. "); perror("malloc"); return 1; }
            cblas_ccopy((int)L,o,0,x1,1);
            cblas_cgemv(CblasColMajor,CblasTrans,(int)L,(int)V,o,X,(int)L,x1,1,z,Y,1);
            cblas_csscal((int)V,ni,Y,1);
            free(x1);
        }
        else if (G==1)
        {
            float *x1;
            if (!(x1=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in mean_c: problem with malloc. "); perror("malloc"); return 1; }
            cblas_ccopy((int)L,o,0,x1,1);
            cblas_cgemv(CblasRowMajor,CblasTrans,(int)L,(int)V,o,X,(int)V,x1,1,z,Y,1);
            cblas_csscal((int)V,ni,Y,1);
            free(x1);
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X+=2)
                {
                    cblas_cdotu_sub((int)L,X,(int)K,o,0,(_Complex float *)Y);
                    *Y++ *= ni; *Y++ *= ni;
                }
            }
        }
    }

    return 0;
}


int mean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni = 1.0 / L;

    if (N==0) {}
    else if (L==1) { cblas_zcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y++ = 0.0; *Y-- = 0.0;
        for (size_t l=0; l<L; l++) { *Y++ += *X++; *Y-- += *X++; }
        *Y++ *= ni; *Y-- *= ni;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;
        const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};

        if (K==1 && (G==1 || B==1))
        {
            double *x1;
            if (!(x1=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in mean_z: problem with malloc. "); perror("malloc"); return 1; }
            cblas_zcopy((int)L,o,0,x1,1);
            cblas_zgemv(CblasColMajor,CblasTrans,(int)L,(int)V,o,X,(int)L,x1,1,z,Y,1);
            cblas_zdscal((int)V,ni,Y,1);
            free(x1);
        }
        else if (G==1)
        {
            double *x1;
            if (!(x1=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in mean_z: problem with malloc. "); perror("malloc"); return 1; }
            cblas_zcopy((int)L,o,0,x1,1);
            cblas_zgemv(CblasRowMajor,CblasTrans,(int)L,(int)V,o,X,(int)V,x1,1,z,Y,1);
            cblas_zdscal((int)V,ni,Y,1);
            free(x1);
        }
        else
        {
            for (size_t g=0; g<G; g++, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; b++, X+=2)
                {
                    cblas_zdotu_sub((int)L,X,(int)K,o,0,(_Complex double *)Y);
                    *Y++ *= ni; *Y++ *= ni;
                }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
