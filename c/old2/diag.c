//Gets kth diagonal of input X as column vector Y

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int diag_s (float *Y, const float *X, const int R, const int C, const char iscolmajor, const int k);
int diag_d (double *Y, const double *X, const int R, const int C, const char iscolmajor, const int k);
int diag_c (float *Y, const float *X, const int R, const int C, const char iscolmajor, const int k);
int diag_z (double *Y, const double *X, const int R, const int C, const char iscolmajor, const int k);

int diag_inplace_s (float *X, const int R, const int C, const char iscolmajor, const int k);
int diag_inplace_d (double *X, const int R, const int C, const char iscolmajor, const int k);
int diag_inplace_c (float *X, const int R, const int C, const char iscolmajor, const int k);
int diag_inplace_z (double *X, const int R, const int C, const char iscolmajor, const int k);


int diag_s (float *Y, const float *X, const int R, const int C, const char iscolmajor, const int k)
{
    int N;

    //Checks
    if (R<1) { fprintf(stderr,"error in diag_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in diag_s: C (ncols X) must be positive\n"); return 1; }

    //Get N (length of output vec)
    if (k>=0 && C>k) { N = (C-k<R) ? C-k : R; }
    else if (k<=0 && R>-k) { N = (R+k<C) ? R+k : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_scopy(N,&X[k*R],R+1,Y,1); }
        else { cblas_scopy(N,&X[-k],R+1,Y,1); }
    }
    else
    {
        if (k>0) { cblas_scopy(N,&X[k],C+1,Y,1); }
        else { cblas_scopy(N,&X[-k*C],C+1,Y,1); }
    }

    return 0;
}


int diag_d (double *Y, const double *X, const int R, const int C, const char iscolmajor, const int k)
{
    int N;

    //Checks
    if (R<1) { fprintf(stderr,"error in diag_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in diag_d: C (ncols X) must be positive\n"); return 1; }
    
    //Get N (length of output vec)
    if (k>=0 && C>k) { N = (C-k<R) ? C-k : R; }
    else if (k<=0 && R>-k) { N = (R+k<C) ? R+k : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_dcopy(N,&X[k*R],R+1,Y,1); }
        else { cblas_dcopy(N,&X[-k],R+1,Y,1); }
    }
    else
    {
        if (k>0) { cblas_dcopy(N,&X[k],C+1,Y,1); }
        else { cblas_dcopy(N,&X[-k*C],C+1,Y,1); }
    }

    return 0;
}


int diag_c (float *Y, const float *X, const int R, const int C, const char iscolmajor, const int k)
{
    int N;

    //Checks
    if (R<1) { fprintf(stderr,"error in diag_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in diag_c: C (ncols X) must be positive\n"); return 1; }
    
    //Get N (length of output vec)
    if (k>=0 && C>k) { N = (C-k<R) ? C-k : R; }
    else if (k<=0 && R>-k) { N = (R+k<C) ? R+k : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_ccopy(N,&X[2*k*R],R+1,Y,1); }
        else { cblas_ccopy(N,&X[-2*k],R+1,Y,1); }
    }
    else
    {
        if (k>0) { cblas_ccopy(N,&X[2*k],C+1,Y,1); }
        else { cblas_ccopy(N,&X[-2*k*C],C+1,Y,1); }
    }

    return 0;
}


int diag_z (double *Y, const double *X, const int R, const int C, const char iscolmajor, const int k)
{
    int N;

    //Checks
    if (R<1) { fprintf(stderr,"error in diag_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in diag_z: C (ncols X) must be positive\n"); return 1; }

    //Get N (length of output vec)
    if (k>=0 && C>k) { N = (C-k<R) ? C-k : R; }
    else if (k<=0 && R>-k) { N = (R+k<C) ? R+k : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_zcopy(N,&X[2*k*R],R+1,Y,1); }
        else { cblas_zcopy(N,&X[-2*k],R+1,Y,1); }
    }
    else
    {
        if (k>0) { cblas_zcopy(N,&X[2*k],C+1,Y,1); }
        else { cblas_zcopy(N,&X[-2*k*C],C+1,Y,1); }
    }

    return 0;
}


int diag_inplace_s (float *X, const int R, const int C, const char iscolmajor, const int k)
{
    int N;

    //Checks
    if (R<1) { fprintf(stderr,"error in diag_inplace_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in diag_inplace_s: C (ncols X) must be positive\n"); return 1; }

    //Get N (length of output vec)
    if (k>=0 && C>k) { N = (C-k<R) ? C-k : R; }
    else if (k<=0 && R>-k) { N = (R+k<C) ? R+k : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_scopy(N,&X[k*R],R+1,X,1); }
        else { cblas_scopy(N,&X[-k],R+1,X,1); }
    }
    else
    {
        if (k>0) { cblas_scopy(N,&X[k],C+1,X,1); }
        else { cblas_scopy(N,&X[-k*C],C+1,X,1); }
    }

    return 0;
}


int diag_inplace_d (double *X, const int R, const int C, const char iscolmajor, const int k)
{
    int N;

    //Checks
    if (R<1) { fprintf(stderr,"error in diag_inplace_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in diag_inplace_d: C (ncols X) must be positive\n"); return 1; }

    //Get N (length of output vec)
    if (k>=0 && C>k) { N = (C-k<R) ? C-k : R; }
    else if (k<=0 && R>-k) { N = (R+k<C) ? R+k : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_dcopy(N,&X[k*R],R+1,X,1); }
        else { cblas_dcopy(N,&X[-k],R+1,X,1); }
    }
    else
    {
        if (k>0) { cblas_dcopy(N,&X[k],C+1,X,1); }
        else { cblas_dcopy(N,&X[-k*C],C+1,X,1); }
    }

    return 0;
}


int diag_inplace_c (float *X, const int R, const int C, const char iscolmajor, const int k)
{
    int N;

    //Checks
    if (R<1) { fprintf(stderr,"error in diag_inplace_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in diag_inplace_c: C (ncols X) must be positive\n"); return 1; }

    //Get N (length of output vec)
    if (k>=0 && C>k) { N = (C-k<R) ? C-k : R; }
    else if (k<=0 && R>-k) { N = (R+k<C) ? R+k : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_ccopy(N,&X[2*k*R],R+1,X,1); }
        else { cblas_ccopy(N,&X[-2*k],R+1,X,1); }
    }
    else
    {
        if (k>0) { cblas_ccopy(N,&X[2*k],C+1,X,1); }
        else { cblas_ccopy(N,&X[-2*k*C],C+1,X,1); }
    }

    return 0;
}


int diag_inplace_z (double *X, const int R, const int C, const char iscolmajor, const int k)
{
    int N;

    //Checks
    if (R<1) { fprintf(stderr,"error in diag_inplace_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in diag_inplace_z: C (ncols X) must be positive\n"); return 1; }

    //Get N (length of output vec)
    if (k>=0 && C>k) { N = (C-k<R) ? C-k : R; }
    else if (k<=0 && R>-k) { N = (R+k<C) ? R+k : C; }
    else { fprintf(stderr,"k out of range [1-R C-1]\n"); return 1; }

    if (iscolmajor)
    {
        if (k>0) { cblas_zcopy(N,&X[2*k*R],R+1,X,1); }
        else { cblas_zcopy(N,&X[-2*k],R+1,X,1); }
    }
    else
    {
        if (k>0) { cblas_zcopy(N,&X[2*k],C+1,X,1); }
        else { cblas_zcopy(N,&X[-2*k*C],C+1,X,1); }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
