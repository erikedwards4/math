//Gets kth diagonal of input X as column vector Y

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int diag_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int k);
int diag_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int k);
int diag_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int k);
int diag_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int k);

int diag_inplace_s (float *X, const size_t R, const size_t C, const char iscolmajor, const int k);
int diag_inplace_d (double *X, const size_t R, const size_t C, const char iscolmajor, const int k);
int diag_inplace_c (float *X, const size_t R, const size_t C, const char iscolmajor, const int k);
int diag_inplace_z (double *X, const size_t R, const size_t C, const char iscolmajor, const int k);


int diag_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int k)
{
    size_t L;  //length of output vec

    if (k>=0 && (int)(int)C>k) { L = (C-(size_t)k<R) ? C-(size_t)k : R; }
    else if (k<=0 && (int)(int)R>-k) { L = ((int)R+k<(int)C) ? (size_t)((int)R+k) : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_scopy((int)L,&X[k*(int)R],(int)R+1,Y,1); }
        else { cblas_scopy((int)L,&X[-k],(int)R+1,Y,1); }
    }
    else
    {
        if (k>0) { cblas_scopy((int)L,&X[k],(int)C+1,Y,1); }
        else { cblas_scopy((int)L,&X[-k*(int)C],(int)C+1,Y,1); }
    }

    return 0;
}


int diag_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int k)
{
    size_t L;  //length of output vec

    if (k>=0 && (int)(int)C>k) { L = (C-(size_t)k<R) ? C-(size_t)k : R; }
    else if (k<=0 && (int)(int)R>-k) { L = ((int)R+k<(int)C) ? (size_t)((int)R+k) : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_dcopy((int)L,&X[k*(int)R],(int)R+1,Y,1); }
        else { cblas_dcopy((int)L,&X[-k],(int)R+1,Y,1); }
    }
    else
    {
        if (k>0) { cblas_dcopy((int)L,&X[k],(int)C+1,Y,1); }
        else { cblas_dcopy((int)L,&X[-k*(int)C],(int)C+1,Y,1); }
    }

    return 0;
}


int diag_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int k)
{
    size_t L;  //length of output vec

    if (k>=0 && (int)C>k) { L = (C-(size_t)k<R) ? C-(size_t)k : R; }
    else if (k<=0 && (int)R>-k) { L = ((int)R+k<(int)C) ? (size_t)((int)R+k) : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_ccopy((int)L,&X[2*(size_t)k*R],(int)R+1,Y,1); }
        else { cblas_ccopy((int)L,&X[(size_t)(-2*k)],(int)R+1,Y,1); }
    }
    else
    {
        if (k>0) { cblas_ccopy((int)L,&X[2*k],(int)C+1,Y,1); }
        else { cblas_ccopy((int)L,&X[(size_t)(-2*k)*C],(int)C+1,Y,1); }
    }

    return 0;
}


int diag_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int k)
{
    size_t L;  //length of output vec

    if (k>=0 && (int)C>k) { L = (C-(size_t)k<R) ? C-(size_t)k : R; }
    else if (k<=0 && (int)R>-k) { L = ((int)R+k<(int)C) ? (size_t)((int)R+k) : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_zcopy((int)L,&X[2*(size_t)k*R],(int)R+1,Y,1); }
        else { cblas_zcopy((int)L,&X[(size_t)(-2*k)],(int)R+1,Y,1); }
    }
    else
    {
        if (k>0) { cblas_zcopy((int)L,&X[2*k],(int)C+1,Y,1); }
        else { cblas_zcopy((int)L,&X[(size_t)(-2*k)*C],(int)C+1,Y,1); }
    }

    return 0;
}


int diag_inplace_s (float *X, const size_t R, const size_t C, const char iscolmajor, const int k)
{
    size_t L;  //length of output vec

    if (k>=0 && (int)C>k) { L = (C-(size_t)k<R) ? C-(size_t)k : R; }
    else if (k<=0 && (int)R>-k) { L = ((int)R+k<(int)C) ? (size_t)((int)R+k) : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_scopy((int)L,&X[(size_t)k*R],(int)R+1,X,1); }
        else { cblas_scopy((int)L,&X[-k],(int)R+1,X,1); }
    }
    else
    {
        if (k>0) { cblas_scopy((int)L,&X[k],(int)C+1,X,1); }
        else { cblas_scopy((int)L,&X[(size_t)(-k)*C],(int)C+1,X,1); }
    }

    return 0;
}


int diag_inplace_d (double *X, const size_t R, const size_t C, const char iscolmajor, const int k)
{
    size_t L;  //length of output vec

    if (k>=0 && (int)C>k) { L = (C-(size_t)k<R) ? C-(size_t)k : R; }
    else if (k<=0 && (int)R>-k) { L = ((int)R+k<(int)C) ? (size_t)((int)R+k) : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_dcopy((int)L,&X[(size_t)k*R],(int)R+1,X,1); }
        else { cblas_dcopy((int)L,&X[-k],(int)R+1,X,1); }
    }
    else
    {
        if (k>0) { cblas_dcopy((int)L,&X[k],(int)C+1,X,1); }
        else { cblas_dcopy((int)L,&X[(size_t)(-k)*C],(int)C+1,X,1); }
    }

    return 0;
}


int diag_inplace_c (float *X, const size_t R, const size_t C, const char iscolmajor, const int k)
{
    size_t L;  //length of output vec

    if (k>=0 && (int)C>k) { L = (C-(size_t)k<R) ? C-(size_t)k : R; }
    else if (k<=0 && (int)R>-k) { L = ((int)R+k<(int)C) ? (size_t)((int)R+k) : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_ccopy((int)L,&X[2*(size_t)k*R],(int)R+1,X,1); }
        else { cblas_ccopy((int)L,&X[(size_t)(-2*k)],(int)R+1,X,1); }
    }
    else
    {
        if (k>0) { cblas_ccopy((int)L,&X[2*k],(int)C+1,X,1); }
        else { cblas_ccopy((int)L,&X[(size_t)(-2*k)*C],(int)C+1,X,1); }
    }

    return 0;
}


int diag_inplace_z (double *X, const size_t R, const size_t C, const char iscolmajor, const int k)
{
    size_t L;  //length of output vec

    if (k>=0 && (int)C>k) { L = (C-(size_t)k<R) ? C-(size_t)k : R; }
    else if (k<=0 && (int)R>-k) { L = ((int)R+k<(int)C) ? (size_t)((int)R+k) : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_zcopy((int)L,&X[2*(size_t)k*R],(int)R+1,X,1); }
        else { cblas_zcopy((int)L,&X[(size_t)(-2*k)],(int)R+1,X,1); }
    }
    else
    {
        if (k>0) { cblas_zcopy((int)L,&X[2*k],(int)C+1,X,1); }
        else { cblas_zcopy((int)L,&X[(size_t)(-2*k)*C],(int)C+1,X,1); }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
