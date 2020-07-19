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

int toeplitz1_s (float *Y, const float *X, const size_t L);
int toeplitz1_d (double *Y, const double *X, const size_t L);
int toeplitz1_c (float *Y, const float *X, const size_t L);
int toeplitz1_z (double *Y, const double *X, const size_t L);

int toeplitz2_s (float *Y, const float *X1, const float *X2, const size_t L1, const size_t L2, const char iscolmajor);
int toeplitz2_d (double *Y, const double *X1, const double *X2, const size_t L1, const size_t L2, const char iscolmajor);
int toeplitz2_c (float *Y, const float *X1, const float *X2, const size_t L1, const size_t L2, const char iscolmajor);
int toeplitz2_z (double *Y, const double *X1, const double *X2, const size_t L1, const size_t L2, const char iscolmajor);



int toeplitz1_s (float *Y, const float *X, const size_t L)
{
    for (size_t l=0; l<L; l++, X++)
    {
        cblas_scopy((int)(L-l),X,0,&Y[l],(int)L+1);
        cblas_scopy((int)(L-l),X,0,&Y[l*L],(int)L+1);
    }

    return 0;
}


int toeplitz1_d (double *Y, const double *X, const size_t L)
{
    for (size_t l=0; l<L; l++, X++)
    {
        cblas_dcopy((int)(L-l),X,0,&Y[l],(int)L+1);
        cblas_dcopy((int)(L-l),X,0,&Y[l*L],(int)L+1);
    }

    return 0;
}


int toeplitz1_c (float *Y, const float *X, const size_t L)
{
    for (size_t l=0; l<L; l++, X+=2)
    {
        cblas_ccopy((int)(L-l),X,0,&Y[2*l],(int)L+1);
        cblas_ccopy((int)(L-l),X,0,&Y[2*l*L],(int)L+1);
    }

    return 0;
}


int toeplitz1_z (double *Y, const double *X, const size_t L)
{
    for (size_t l=0; l<L; l++, X+=2)
    {
        cblas_zcopy((int)(L-l),X,0,&Y[2*l],(int)L+1);
        cblas_zcopy((int)(L-l),X,0,&Y[2*l*L],(int)L+1);
    }

    return 0;
}


int toeplitz2_s (float *Y, const float *X1, const float *X2, const size_t L1, const size_t L2, const char iscolmajor)
{
    if (iscolmajor)
    {
        for (size_t l=0; l<L1; l++, X1++) { cblas_scopy((int)(L1-l),X1,0,&Y[l],(int)L1+1); }
        for (size_t l=0; l<L2; l++, X2++) { cblas_scopy((int)(L2-l),X2,0,&Y[l*L1],(int)L1+1); }
    }
    else
    {
        for (size_t l=0; l<L1; l++, X1++) { cblas_scopy((int)(L1-l),X1,0,&Y[l*L2],(int)L2+1); }
        for (size_t l=0; l<L2; l++, X2++) { cblas_scopy((int)(L2-l),X2,0,&Y[l],(int)L2+1); }   
    }

    return 0;
}


int toeplitz2_d (double *Y, const double *X1, const double *X2, const size_t L1, const size_t L2, const char iscolmajor)
{
    if (iscolmajor)
    {
        for (size_t l=0; l<L1; l++, X1++) { cblas_dcopy((int)(L1-l),X1,0,&Y[l],(int)L1+1); }
        for (size_t l=0; l<L2; l++, X2++) { cblas_dcopy((int)(L2-l),X2,0,&Y[l*L1],(int)L1+1); }
    }
    else
    {
        for (size_t l=0; l<L1; l++, X1++) { cblas_dcopy((int)(L1-l),X1,0,&Y[l*L2],(int)L2+1); }
        for (size_t l=0; l<L2; l++, X2++) { cblas_dcopy((int)(L2-l),X2,0,&Y[l],(int)L2+1); }
    }

    return 0;
}


int toeplitz2_c (float *Y, const float *X1, const float *X2, const size_t L1, const size_t L2, const char iscolmajor)
{
    if (iscolmajor)
    {
        for (size_t l=0; l<L1; l++, X1+=2) { cblas_ccopy((int)(L1-l),X1,0,&Y[2*l],(int)L1+1); }
        for (size_t l=0; l<L2; l++, X2+=2) { cblas_ccopy((int)(L2-l),X2,0,&Y[2*l*L1],(int)L1+1); }
    }
    else
    {
        for (size_t l=0; l<L1; l++, X1+=2) { cblas_ccopy((int)(L1-l),X1,0,&Y[2*l*L2],(int)L2+1); }
        for (size_t l=0; l<L2; l++, X2+=2) { cblas_ccopy((int)(L2-l),X2,0,&Y[2*l],(int)L2+1); }
    }

    return 0;
}


int toeplitz2_z (double *Y, const double *X1, const double *X2, const size_t L1, const size_t L2, const char iscolmajor)
{
    if (iscolmajor)
    {
        for (size_t l=0; l<L1; l++, X1+=2) { cblas_zcopy((int)(L1-l),X1,0,&Y[2*l],(int)L1+1); }
        for (size_t l=0; l<L2; l++, X2+=2) { cblas_zcopy((int)(L2-l),X2,0,&Y[2*l*L1],(int)L1+1); }
    }
    else
    {
        for (size_t l=0; l<L1; l++, X1+=2) { cblas_zcopy((int)(L1-l),X1,0,&Y[2*l*L2],(int)L2+1); }
        for (size_t l=0; l<L2; l++, X2+=2) { cblas_zcopy((int)(L2-l),X2,0,&Y[2*l],(int)L2+1); }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
