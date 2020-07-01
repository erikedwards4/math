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

int toeplitz1_s (float *Y, const float *X, const int N);
int toeplitz1_d (double *Y, const double *X, const int N);
int toeplitz1_c (float *Y, const float *X, const int N);
int toeplitz1_z (double *Y, const double *X, const int N);

int toeplitz2_s (float *Y, const float *X1, const float *X2, const int N1, const int N2, const char iscolmajor);
int toeplitz2_d (double *Y, const double *X1, const double *X2, const int N1, const int N2, const char iscolmajor);
int toeplitz2_c (float *Y, const float *X1, const float *X2, const int N1, const int N2, const char iscolmajor);
int toeplitz2_z (double *Y, const double *X1, const double *X2, const int N1, const int N2, const char iscolmajor);



int toeplitz1_s (float *Y, const float *X, const int N)
{
    if (N<1) { fprintf(stderr,"error in toeplitz1_s: N (length X1) must be positive\n"); return 1; }

    for (int k=0; k<N; k++)
    {
        cblas_scopy(N-k,&X[k],0,&Y[k],N+1);
        cblas_scopy(N-k,&X[k],0,&Y[k*N],N+1);
    }

    return 0;
}


int toeplitz1_d (double *Y, const double *X, const int N)
{
    if (N<1) { fprintf(stderr,"error in toeplitz1_d: N (length X1) must be positive\n"); return 1; }

    for (int k=0; k<N; k++)
    {
        cblas_dcopy(N-k,&X[k],0,&Y[k],N+1);
        cblas_dcopy(N-k,&X[k],0,&Y[k*N],N+1);
    }

    return 0;
}


int toeplitz1_c (float *Y, const float *X, const int N)
{
    if (N<1) { fprintf(stderr,"error in toeplitz1_c: N (length X1) must be positive\n"); return 1; }

    for (int k=0; k<N; k++)
    {
        cblas_ccopy(N-k,&X[2*k],0,&Y[2*k],N+1);
        cblas_ccopy(N-k,&X[2*k],0,&Y[2*k*N],N+1);
    }

    return 0;
}


int toeplitz1_z (double *Y, const double *X, const int N)
{
    if (N<1) { fprintf(stderr,"error in toeplitz1_z: N (length X1) must be positive\n"); return 1; }

    for (int k=0; k<N; k++)
    {
        cblas_zcopy(N-k,&X[2*k],0,&Y[2*k],N+1);
        cblas_zcopy(N-k,&X[2*k],0,&Y[2*k*N],N+1);
    }

    return 0;
}


int toeplitz2_s (float *Y, const float *X1, const float *X2, const int N1, const int N2, const char iscolmajor)
{
    if (N1<1) { fprintf(stderr,"error in toeplitz2_s: N1 (length X1) must be positive\n"); return 1; }
    if (N2<1) { fprintf(stderr,"error in toeplitz2_s: N2 (length X2) must be positive\n"); return 1; }

    int k;

    if (iscolmajor)
    {
        for (k=0; k<N1; k++) { cblas_scopy(N1-k,&X1[k],0,&Y[k],N1+1); }
        for (k=0; k<N2; k++) { cblas_scopy(N2-k,&X2[k],0,&Y[k*N1],N1+1); }
    }
    else
    {
        for (k=0; k<N1; k++) { cblas_scopy(N1-k,&X1[k],0,&Y[k*N2],N2+1); }
        for (k=0; k<N2; k++) { cblas_scopy(N2-k,&X2[k],0,&Y[k],N2+1); }   
    }

    return 0;
}


int toeplitz2_d (double *Y, const double *X1, const double *X2, const int N1, const int N2, const char iscolmajor)
{
    if (N1<1) { fprintf(stderr,"error in toeplitz2_d: N1 (length X1) must be positive\n"); return 1; }
    if (N2<1) { fprintf(stderr,"error in toeplitz2_d: N2 (length X2) must be positive\n"); return 1; }

    int k;

    if (iscolmajor)
    {
        for (k=0; k<N1; k++) { cblas_dcopy(N1-k,&X1[k],0,&Y[k],N1+1); }
        for (k=0; k<N2; k++) { cblas_dcopy(N2-k,&X2[k],0,&Y[k*N1],N1+1); }
    }
    else
    {
        for (k=0; k<N1; k++) { cblas_dcopy(N1-k,&X1[k],0,&Y[k*N2],N2+1); }
        for (k=0; k<N2; k++) { cblas_dcopy(N2-k,&X2[k],0,&Y[k],N2+1); }
    }

    return 0;
}


int toeplitz2_c (float *Y, const float *X1, const float *X2, const int N1, const int N2, const char iscolmajor)
{
    if (N1<1) { fprintf(stderr,"error in toeplitz2_c: N1 (length X1) must be positive\n"); return 1; }
    if (N2<1) { fprintf(stderr,"error in toeplitz2_c: N2 (length X2) must be positive\n"); return 1; }

    int k;

    if (iscolmajor)
    {
        for (k=0; k<N1; k++) { cblas_ccopy(N1-k,&X1[2*k],0,&Y[2*k],N1+1); }
        for (k=0; k<N2; k++) { cblas_ccopy(N2-k,&X2[2*k],0,&Y[2*k*N1],N1+1); }
    }
    else
    {
        for (k=0; k<N1; k++) { cblas_ccopy(N1-k,&X1[2*k],0,&Y[2*k*N2],N2+1); }
        for (k=0; k<N2; k++) { cblas_ccopy(N2-k,&X2[2*k],0,&Y[2*k],N2+1); }
    }

    return 0;
}


int toeplitz2_z (double *Y, const double *X1, const double *X2, const int N1, const int N2, const char iscolmajor)
{
    if (N1<1) { fprintf(stderr,"error in toeplitz2_z: N1 (length X1) must be positive\n"); return 1; }
    if (N2<1) { fprintf(stderr,"error in toeplitz2_z: N2 (length X2) must be positive\n"); return 1; }

    int k;

    if (iscolmajor)
    {
        for (k=0; k<N1; k++) { cblas_zcopy(N1-k,&X1[2*k],0,&Y[2*k],N1+1); }
        for (k=0; k<N2; k++) { cblas_zcopy(N2-k,&X2[2*k],0,&Y[2*k*N1],N1+1); }
    }
    else
    {
        for (k=0; k<N1; k++) { cblas_zcopy(N1-k,&X1[2*k],0,&Y[2*k*N2],N2+1); }
        for (k=0; k<N2; k++) { cblas_zcopy(N2-k,&X2[2*k],0,&Y[2*k],N2+1); }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
