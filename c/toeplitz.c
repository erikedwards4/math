//Makes Toeplitz matrix Y from vec X1 (1st col of Y), and from vec X2 (1st row of Y).
//
//This has versions toeplitz1 and toeplitz2 for 1 input (X1)
//or 2 inputs (X1,X2). In toeplitz1, it is implied that X2==X1.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int toeplitz1_s (float *Y, const float *X, const size_t N);
int toeplitz1_d (double *Y, const double *X, const size_t N);
int toeplitz1_c (float *Y, const float *X, const size_t N);
int toeplitz1_z (double *Y, const double *X, const size_t N);

int toeplitz2_s (float *Y, const float *X1, const float *X2, const size_t N1, const size_t N2, const char iscolmajor);
int toeplitz2_d (double *Y, const double *X1, const double *X2, const size_t N1, const size_t N2, const char iscolmajor);
int toeplitz2_c (float *Y, const float *X1, const float *X2, const size_t N1, const size_t N2, const char iscolmajor);
int toeplitz2_z (double *Y, const double *X1, const double *X2, const size_t N1, const size_t N2, const char iscolmajor);



int toeplitz1_s (float *Y, const float *X, const size_t N)
{
    for (size_t k=0; k<N; k++)
    {
        cblas_scopy((int)(N-k),&X[k],0,&Y[k],(int)N+1);
        cblas_scopy((int)(N-k),&X[k],0,&Y[k*N],(int)N+1);
    }

    return 0;
}


int toeplitz1_d (double *Y, const double *X, const size_t N)
{
    for (size_t k=0; k<N; k++)
    {
        cblas_dcopy((int)(N-k),&X[k],0,&Y[k],(int)N+1);
        cblas_dcopy((int)(N-k),&X[k],0,&Y[k*N],(int)N+1);
    }

    return 0;
}


int toeplitz1_c (float *Y, const float *X, const size_t N)
{
    for (size_t k=0; k<N; k++)
    {
        cblas_ccopy((int)(N-k),&X[2*k],0,&Y[2*k],(int)N+1);
        cblas_ccopy((int)(N-k),&X[2*k],0,&Y[2*k*N],(int)N+1);
    }

    return 0;
}


int toeplitz1_z (double *Y, const double *X, const size_t N)
{
    for (size_t k=0; k<N; k++)
    {
        cblas_zcopy((int)(N-k),&X[2*k],0,&Y[2*k],(int)N+1);
        cblas_zcopy((int)(N-k),&X[2*k],0,&Y[2*k*N],(int)N+1);
    }

    return 0;
}


int toeplitz2_s (float *Y, const float *X1, const float *X2, const size_t N1, const size_t N2, const char iscolmajor)
{
    if (iscolmajor)
    {
        for (size_t k=0; k<N1; k++) { cblas_scopy((int)(N1-k),&X1[k],0,&Y[k],(int)N1+1); }
        for (size_t k=0; k<N2; k++) { cblas_scopy((int)(N2-k),&X2[k],0,&Y[k*N1],(int)N1+1); }
    }
    else
    {
        for (size_t k=0; k<N1; k++) { cblas_scopy((int)(N1-k),&X1[k],0,&Y[k*N2],(int)N2+1); }
        for (size_t k=0; k<N2; k++) { cblas_scopy((int)(N2-k),&X2[k],0,&Y[k],(int)N2+1); }   
    }

    return 0;
}


int toeplitz2_d (double *Y, const double *X1, const double *X2, const size_t N1, const size_t N2, const char iscolmajor)
{
    if (iscolmajor)
    {
        for (size_t k=0; k<N1; k++) { cblas_dcopy((int)(N1-k),&X1[k],0,&Y[k],(int)N1+1); }
        for (size_t k=0; k<N2; k++) { cblas_dcopy((int)(N2-k),&X2[k],0,&Y[k*N1],(int)N1+1); }
    }
    else
    {
        for (size_t k=0; k<N1; k++) { cblas_dcopy((int)(N1-k),&X1[k],0,&Y[k*N2],(int)N2+1); }
        for (size_t k=0; k<N2; k++) { cblas_dcopy((int)(N2-k),&X2[k],0,&Y[k],(int)N2+1); }
    }

    return 0;
}


int toeplitz2_c (float *Y, const float *X1, const float *X2, const size_t N1, const size_t N2, const char iscolmajor)
{
    if (iscolmajor)
    {
        for (size_t k=0; k<N1; k++) { cblas_ccopy((int)(N1-k),&X1[2*k],0,&Y[2*k],(int)N1+1); }
        for (size_t k=0; k<N2; k++) { cblas_ccopy((int)(N2-k),&X2[2*k],0,&Y[2*k*N1],(int)N1+1); }
    }
    else
    {
        for (size_t k=0; k<N1; k++) { cblas_ccopy((int)(N1-k),&X1[2*k],0,&Y[2*k*N2],(int)N2+1); }
        for (size_t k=0; k<N2; k++) { cblas_ccopy((int)(N2-k),&X2[2*k],0,&Y[2*k],(int)N2+1); }
    }

    return 0;
}


int toeplitz2_z (double *Y, const double *X1, const double *X2, const size_t N1, const size_t N2, const char iscolmajor)
{
    if (iscolmajor)
    {
        for (size_t k=0; k<N1; k++) { cblas_zcopy((int)(N1-k),&X1[2*k],0,&Y[2*k],(int)N1+1); }
        for (size_t k=0; k<N2; k++) { cblas_zcopy((int)(N2-k),&X2[2*k],0,&Y[2*k*N1],(int)N1+1); }
    }
    else
    {
        for (size_t k=0; k<N1; k++) { cblas_zcopy((int)(N1-k),&X1[2*k],0,&Y[2*k*N2],(int)N2+1); }
        for (size_t k=0; k<N2; k++) { cblas_zcopy((int)(N2-k),&X2[2*k],0,&Y[2*k],(int)N2+1); }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
