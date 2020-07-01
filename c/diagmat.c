//Puts vector X on kth diagonal of matrix Y.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int diagmat_s (float *Y, const float *X, const int R, const int C, const char iscolmajor, const int k);
int diagmat_d (double *Y, const double *X, const int R, const int C, const char iscolmajor, const int k);
int diagmat_c (float *Y, const float *X, const int R, const int C, const char iscolmajor, const int k);
int diagmat_z (double *Y, const double *X, const int R, const int C, const char iscolmajor, const int k);



int diagmat_s (float *Y, const float *X, const int R, const int C, const char iscolmajor, const int k)
{
    if (R<1) { fprintf(stderr,"error in diagmat_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in diagmat_s: C (ncols X) must be positive\n"); return 1; }
    
    const int N = (k>0) ? R : C;
    const float z = 0.0f;

    cblas_scopy(R*C,&z,0,Y,1);

    if (iscolmajor)
    {
        if (k<0) { cblas_scopy(N,X,1,&Y[-k],R+1); }
        else { cblas_scopy(N,X,1,&Y[k*R],R+1); }
    }
    else
    {
        if (k<0) { cblas_scopy(N,X,1,&Y[-k*C],C+1); }
        else { cblas_scopy(N,X,1,&Y[k],C+1); }
    }

    return 0;
}


int diagmat_d (double *Y, const double *X, const int R, const int C, const char iscolmajor, const int k)
{
    if (R<1) { fprintf(stderr,"error in diagmat_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in diagmat_d: C (ncols X) must be positive\n"); return 1; }

    const int N = (k>0) ? R : C;
    const double z = 0.0;

    cblas_dcopy(R*C,&z,0,Y,1);

    if (iscolmajor)
    {
        if (k<0) { cblas_dcopy(N,X,1,&Y[-k],R+1); }
        else { cblas_dcopy(N,X,1,&Y[k*R],R+1); }
    }
    else
    {
        if (k<0) { cblas_dcopy(N,X,1,&Y[-k*C],C+1); }
        else { cblas_dcopy(N,X,1,&Y[k],C+1); }
    }

    return 0;
}


int diagmat_c (float *Y, const float *X, const int R, const int C, const char iscolmajor, const int k)
{
    if (R<1) { fprintf(stderr,"error in diagmat_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in diagmat_c: C (ncols X) must be positive\n"); return 1; }

    const int N = (k>0) ? R : C;
    const float z[2] = {0.0f,0.0f};

    cblas_ccopy(R*C,z,0,Y,1);

    if (iscolmajor)
    {
        if (k<0) { cblas_ccopy(N,X,1,&Y[-2*k],R+1); }
        else { cblas_ccopy(N,X,1,&Y[2*k*R],R+1); }
    }
    else
    {
        if (k<0) { cblas_ccopy(N,X,1,&Y[-2*k*C],C+1); }
        else { cblas_ccopy(N,X,1,&Y[2*k],C+1); }
    }

    return 0;
}


int diagmat_z (double *Y, const double *X, const int R, const int C, const char iscolmajor, const int k)
{
    if (R<1) { fprintf(stderr,"error in diagmat_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in diagmat_z: C (ncols X) must be positive\n"); return 1; }

    const int N = (k>0) ? R : C;
    const double z[2] = {0.0,0.0};

    cblas_zcopy(R*C,z,0,Y,1);

    if (iscolmajor)
    {
        if (k<0) { cblas_zcopy(N,X,1,&Y[-2*k],R+1); }
        else { cblas_zcopy(N,X,1,&Y[2*k*R],R+1); }
    }
    else
    {
        if (k<0) { cblas_zcopy(N,X,1,&Y[-2*k*C],C+1); }
        else { cblas_zcopy(N,X,1,&Y[2*k],C+1); }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
