//Gets one column of input X

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int col_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c);
int col_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c);
int col_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c);
int col_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c);

int col_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c);
int col_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c);
int col_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c);
int col_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c);


int col_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c)
{
    const int HS = H*S;
    int r, s, h, m, n = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in col_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in col_s: C (ncols X) must be positive\n"); return 1; }
    if (C<=c) { fprintf(stderr,"error in col_s: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (iscolmajor)
    {
        m = c*R;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                cblas_scopy(R,&X[m],1,&Y[n],1);
                m += R*C; n += R;
            }
        }
    }
    else if (HS==1)
    {
        cblas_scopy(R,&X[c],C,Y,1);
    }
    else
    {
        for (r=0; r<R; r++)
        {
            cblas_scopy(HS,&X[c*HS+r*HS*C],1,&Y[r*HS],1);
        }
    }

    return 0;
}


int col_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c)
{
    const int HS = H*S;
    int r, s, h, m, n = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in col_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in col_d: C (ncols X) must be positive\n"); return 1; }
    if (C<=c) { fprintf(stderr,"error in col_d: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (iscolmajor)
    {
        m = c*R;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                cblas_dcopy(R,&X[m],1,&Y[n],1);
                m += R*C; n += R;
            }
        }
    }
    else if (HS==1)
    {
        cblas_dcopy(R,&X[c],C,Y,1);
    }
    else
    {
        for (r=0; r<R; r++)
        {
            cblas_dcopy(HS,&X[c*HS+r*HS*C],1,&Y[r*HS],1);
        }
    }

    return 0;
}


int col_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c)
{
    const int HS = H*S;
    int r, s, h, m, n = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in col_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in col_c: C (ncols X) must be positive\n"); return 1; }
    if (C<=c) { fprintf(stderr,"error in col_c: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (iscolmajor)
    {
        m = 2*c*R;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                cblas_ccopy(R,&X[m],1,&Y[n],1);
                m += 2*R*C; n += 2*R;
            }
        }
    }
    else if (HS==1)
    {
        cblas_ccopy(R,&X[2*c],C,Y,1);
    }
    else
    {
        for (r=0; r<R; r++)
        {
            cblas_ccopy(HS,&X[2*(c*HS+r*HS*C)],1,&Y[2*r*HS],1);
        }
    }

    return 0;
}


int col_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c)
{
    const int HS = H*S;
    int r, s, h, m, n = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in col_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in col_z: C (ncols X) must be positive\n"); return 1; }
    if (C<=c) { fprintf(stderr,"error in col_z: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (iscolmajor)
    {
        m = 2*c*R;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                cblas_zcopy(R,&X[m],1,&Y[n],1);
                m += 2*R*C; n += 2*R;
            }
        }
    }
    else if (HS==1)
    {
        cblas_zcopy(R,&X[2*c],C,Y,1);
    }
    else
    {
        for (r=0; r<R; r++)
        {
            cblas_zcopy(HS,&X[2*(c*HS+r*HS*C)],1,&Y[2*r*HS],1);
        }
    }

    return 0;
}


int col_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c)
{
    int s, h, m, n = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in col_inplace_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in col_inplace_s: C (ncols X) must be positive\n"); return 1; }
    if (C<=c) { fprintf(stderr,"error in col_inplace_s: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (iscolmajor)
    {
        m = c*R;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                if (m!=n) { cblas_scopy(R,&X[m],1,&X[n],1); }
                m += R*C; n += R;
            }
        }
    }
    else if (c>0)
    {
        const int HS = H*S;
        if (HS==1) { cblas_scopy(R,&X[c],C,X,1); }
        else
        {
            for (int r=0; r<R; r++)
            {
                cblas_scopy(HS,&X[c*HS+r*HS*C],1,&X[r*HS],1);
            }
        }
    }

    return 0;
}


int col_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c)
{
    int s, h, m, n = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in col_inplace_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in col_inplace_d: C (ncols X) must be positive\n"); return 1; }
    if (C<=c) { fprintf(stderr,"error in col_inplace_d: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (iscolmajor)
    {
        m = c*R;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                if (m!=n) { cblas_dcopy(R,&X[m],1,&X[n],1); }
                m += R*C; n += R;
            }
        }
    }
    else if (c>0)
    {
        const int HS = H*S;
        if (HS==1) { cblas_dcopy(R,&X[c],C,X,1); }
        else
        {
            for (int r=0; r<R; r++)
            {
                cblas_dcopy(HS,&X[c*HS+r*HS*C],1,&X[r*HS],1);
            }
        }
    }

    return 0;
}


int col_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c)
{
    int s, h, m, n = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in col_inplace_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in col_inplace_c: C (ncols X) must be positive\n"); return 1; }
    if (C<=c) { fprintf(stderr,"error in col_inplace_c: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (iscolmajor)
    {
        m = 2*c*R;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                if (m!=n) { cblas_ccopy(R,&X[m],1,&X[n],1); }
                m += 2*R*C; n += 2*R;
            }
        }
    }
    else if (c>0)
    {
        const int HS = H*S;
        if (HS==1) { cblas_ccopy(R,&X[2*c],C,X,1); }
        else
        {
            for (int r=0; r<R; r++)
            {
                cblas_ccopy(HS,&X[2*(c*HS+r*HS*C)],1,&X[2*r*HS],1);
            }
        }
    }

    return 0;
}


int col_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int c)
{
    int s, h, m, n = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in col_inplace_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in col_inplace_z: C (ncols X) must be positive\n"); return 1; }
    if (C<=c) { fprintf(stderr,"error in col_inplace_z: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (iscolmajor)
    {
        m = 2*c*R;
        for (h=0; h<H; h++)
        {
            for (s=0; s<S; s++)
            {
                if (m!=n) { cblas_zcopy(R,&X[m],1,&X[n],1); }
                m += 2*R*C; n += 2*R;
            }
        }
    }
    else if (c>0)
    {
        const int HS = H*S;
        if (HS==1) { cblas_zcopy(R,&X[2*c],C,X,1); }
        else
        {
            for (int r=0; r<R; r++)
            {
                cblas_zcopy(HS,&X[2*(c*HS+r*HS*C)],1,&X[2*r*HS],1);
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
