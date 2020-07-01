//Gets one row of input X

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int row_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r);
int row_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r);
int row_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r);
int row_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r);

int row_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r);
int row_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r);
int row_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r);
int row_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r);


int row_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r)
{
    if (R<1) { fprintf(stderr,"error in row_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in row_s: C (ncols X) must be positive\n"); return 1; }
    if (R<=r) { fprintf(stderr,"error in row_s: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        int n = 0, m = r;
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                cblas_scopy(C,&X[m],R,&Y[n],1);
                m += R*C; n += R;
            }
        }
    }
    else if (r>0)
    {
        cblas_scopy(H*S*C,&X[r*H*S*C],1,Y,1);
    }

    return 0;
}


int row_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r)
{
    if (R<1) { fprintf(stderr,"error in row_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in row_d: C (ncols X) must be positive\n"); return 1; }
    if (R<=r) { fprintf(stderr,"error in row_d: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        int n = 0, m = r;
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                cblas_dcopy(C,&X[m],R,&Y[n],1);
                m += R*C; n += R;
            }
        }
    }
    else if (r>0)
    {
        cblas_dcopy(H*S*C,&X[r*H*S*C],1,Y,1);
    }

    return 0;
}


int row_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r)
{
    if (R<1) { fprintf(stderr,"error in row_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in row_c: C (ncols X) must be positive\n"); return 1; }
    if (R<=r) { fprintf(stderr,"error in row_c: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        int n = 0, m = 2*r;
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                cblas_ccopy(C,&X[m],R,&Y[n],1);
                m += 2*R*C; n += 2*R;
            }
        }
    }
    else if (r>0)
    {
        cblas_ccopy(H*S*C,&X[2*r*H*S*C],1,Y,1);
    }

    return 0;
}


int row_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r)
{
    if (R<1) { fprintf(stderr,"error in row_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in row_z: C (ncols X) must be positive\n"); return 1; }
    if (R<=r) { fprintf(stderr,"error in row_z: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        int n = 0, m = 2*r;
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                cblas_zcopy(C,&X[m],R,&Y[n],1);
                m += 2*R*C; n += 2*R;
            }
        }
    }
    else if (r>0)
    {
        cblas_zcopy(H*S*C,&X[2*r*H*S*C],1,Y,1);
    }

    return 0;
}


int row_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r)
{
    if (R<1) { fprintf(stderr,"error in row_inplace_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in row_inplace_s: C (ncols X) must be positive\n"); return 1; }
    if (R<=r) { fprintf(stderr,"error in row_inplace_s: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        int n = 0, m = r;
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                cblas_scopy(C,&X[m],R,&X[n],1);
                m += R*C; n += R;
            }
        }
    }
    else if (r>0)
    {
        cblas_scopy(H*S*C,&X[r*H*S*C],1,X,1);
    }

    return 0;
}


int row_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r)
{
    if (R<1) { fprintf(stderr,"error in row_inplace_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in row_inplace_d: C (ncols X) must be positive\n"); return 1; }
    if (R<=r) { fprintf(stderr,"error in row_inplace_d: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        int n = 0, m = r;
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                cblas_dcopy(C,&X[m],R,&X[n],1);
                m += R*C; n += R;
            }
        }
    }
    else if (r>0)
    {
        cblas_dcopy(H*S*C,&X[r*H*S*C],1,X,1);
    }

    return 0;
}


int row_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r)
{
    if (R<1) { fprintf(stderr,"error in row_inplace_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in row_inplace_c: C (ncols X) must be positive\n"); return 1; }
    if (R<=r) { fprintf(stderr,"error in row_inplace_c: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        int n = 0, m = 2*r;
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                cblas_ccopy(C,&X[m],R,&X[n],1);
                m += 2*R*C; n += 2*R;
            }
        }
    }
    else if (r>0)
    {
        cblas_ccopy(H*S*C,&X[2*r*H*S*C],1,X,1);
    }

    return 0;
}


int row_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int r)
{
    if (R<1) { fprintf(stderr,"error in row_inplace_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in row_inplace_z: C (ncols X) must be positive\n"); return 1; }
    if (R<=r) { fprintf(stderr,"error in row_inplace_z: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        int n = 0, m = 2*r;
        for (int h=0; h<H; h++)
        {
            for (int s=0; s<S; s++)
            {
                cblas_zcopy(C,&X[m],R,&X[n],1);
                m += 2*R*C; n += 2*R;
            }
        }
    }
    else if (r>0)
    {
        cblas_zcopy(H*S*C,&X[2*r*H*S*C],1,X,1);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
