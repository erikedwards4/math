//Puts vector X on kth diagonal of matrix Y.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int diagmat_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int k);
int diagmat_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int k);
int diagmat_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int k);
int diagmat_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int k);



int diagmat_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int k)
{
    const size_t N = (k>0) ? R : C;
    const float z = 0.0f;

    cblas_scopy((int)(R*C),&z,0,Y,1);

    if (iscolmajor)
    {
        if (k<0) { cblas_scopy((int)N,X,1,&Y[-k],(int)R+1); }
        else { cblas_scopy((int)N,X,1,&Y[(size_t)k*R],(int)R+1); }
    }
    else
    {
        if (k<0) { cblas_scopy((int)N,X,1,&Y[(size_t)(-k)*C],(int)C+1); }
        else { cblas_scopy((int)N,X,1,&Y[k],(int)C+1); }
    }

    return 0;
}


int diagmat_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int k)
{
    const size_t N = (k>0) ? R : C;
    const double z = 0.0;

    cblas_dcopy((int)(R*C),&z,0,Y,1);

    if (iscolmajor)
    {
        if (k<0) { cblas_dcopy((int)N,X,1,&Y[-k],(int)R+1); }
        else { cblas_dcopy((int)N,X,1,&Y[(size_t)k*R],(int)R+1); }
    }
    else
    {
        if (k<0) { cblas_dcopy((int)N,X,1,&Y[(size_t)(-k)*C],(int)C+1); }
        else { cblas_dcopy((int)N,X,1,&Y[k],(int)C+1); }
    }

    return 0;
}


int diagmat_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int k)
{
    const size_t N = (k>0) ? R : C;
    const float z[2] = {0.0f,0.0f};

    cblas_ccopy((int)(R*C),z,0,Y,1);

    if (iscolmajor)
    {
        if (k<0) { cblas_ccopy((int)N,X,1,&Y[-2*k],(int)R+1); }
        else { cblas_ccopy((int)N,X,1,&Y[(size_t)(2*k)*R],(int)R+1); }
    }
    else
    {
        if (k<0) { cblas_ccopy((int)N,X,1,&Y[(size_t)(-2*k)*C],(int)C+1); }
        else { cblas_ccopy((int)N,X,1,&Y[2*k],(int)C+1); }
    }

    return 0;
}


int diagmat_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int k)
{
    const size_t N = (k>0) ? R : C;
    const double z[2] = {0.0,0.0};

    cblas_zcopy((int)(R*C),z,0,Y,1);

    if (iscolmajor)
    {
        if (k<0) { cblas_zcopy((int)N,X,1,&Y[-2*k],(int)R+1); }
        else { cblas_zcopy((int)N,X,1,&Y[(size_t)(2*k)*R],(int)R+1); }
    }
    else
    {
        if (k<0) { cblas_zcopy((int)N,X,1,&Y[(size_t)(-2*k)*C],(int)C+1); }
        else { cblas_zcopy((int)N,X,1,&Y[2*k],(int)C+1); }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
