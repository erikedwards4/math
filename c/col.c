//Gets one column of input X

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int col_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c);
int col_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c);
int col_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c);
int col_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c);

int col_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c);
int col_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c);
int col_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c);
int col_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c);


int col_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c)
{
    if (C<=c) { fprintf(stderr,"error in col_s: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    const size_t HS = H*S;

    if (iscolmajor)
    {
        size_t m = c*R, n = 0;
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                cblas_scopy((int)R,&X[m],1,&Y[n],1);
                m += R*C; n += R;
            }
        }
    }
    else if (HS==1)
    {
        cblas_scopy((int)R,&X[c],(int)C,Y,1);
    }
    else
    {
        for (size_t r=0; r<R; r++)
        {
            cblas_scopy((int)HS,&X[c*HS+r*HS*C],1,&Y[r*HS],1);
        }
    }

    return 0;
}


int col_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c)
{
    if (C<=c) { fprintf(stderr,"error in col_d: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    const size_t HS = H*S;

    if (iscolmajor)
    {
        size_t m = c*R, n = 0;
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                cblas_dcopy((int)R,&X[m],1,&Y[n],1);
                m += R*C; n += R;
            }
        }
    }
    else if (HS==1)
    {
        cblas_dcopy((int)R,&X[c],(int)C,Y,1);
    }
    else
    {
        for (size_t r=0; r<R; r++)
        {
            cblas_dcopy((int)HS,&X[c*HS+r*HS*C],1,&Y[r*HS],1);
        }
    }

    return 0;
}


int col_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c)
{
    if (C<=c) { fprintf(stderr,"error in col_c: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    const size_t HS = H*S;

    if (iscolmajor)
    {
        size_t m = 2*c*R, n = 0;
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                cblas_ccopy((int)R,&X[m],1,&Y[n],1);
                m += 2*R*C; n += 2*R;
            }
        }
    }
    else if (HS==1)
    {
        cblas_ccopy((int)R,&X[2*c],(int)C,Y,1);
    }
    else
    {
        for (size_t r=0; r<R; r++)
        {
            cblas_ccopy((int)HS,&X[2*(c*HS+r*HS*C)],1,&Y[2*r*HS],1);
        }
    }

    return 0;
}


int col_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c)
{
    if (C<=c) { fprintf(stderr,"error in col_z: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    const size_t HS = H*S;

    if (iscolmajor)
    {
        size_t m = 2*c*R, n = 0;
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                cblas_zcopy((int)R,&X[m],1,&Y[n],1);
                m += 2*R*C; n += 2*R;
            }
        }
    }
    else if (HS==1)
    {
        cblas_zcopy((int)R,&X[2*c],(int)C,Y,1);
    }
    else
    {
        for (size_t r=0; r<R; r++)
        {
            cblas_zcopy((int)HS,&X[2*(c*HS+r*HS*C)],1,&Y[2*r*HS],1);
        }
    }

    return 0;
}


int col_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c)
{
    if (C<=c) { fprintf(stderr,"error in col_inplace_s: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (iscolmajor)
    {
        size_t m = c*R, n = 0;
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                if (m!=n) { cblas_scopy((int)R,&X[m],1,&X[n],1); }
                m += R*C; n += R;
            }
        }
    }
    else if (c>0)
    {
        const size_t HS = H*S;
        if (HS==1) { cblas_scopy((int)R,&X[c],(int)C,X,1); }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_scopy((int)HS,&X[c*HS+r*HS*C],1,&X[r*HS],1);
            }
        }
    }

    return 0;
}


int col_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c)
{
    if (C<=c) { fprintf(stderr,"error in col_inplace_d: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (iscolmajor)
    {
        size_t m = c*R, n = 0;
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                if (m!=n) { cblas_dcopy((int)R,&X[m],1,&X[n],1); }
                m += R*C; n += R;
            }
        }
    }
    else if (c>0)
    {
        const size_t HS = H*S;
        if (HS==1) { cblas_dcopy((int)R,&X[c],(int)C,X,1); }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_dcopy((int)HS,&X[c*HS+r*HS*C],1,&X[r*HS],1);
            }
        }
    }

    return 0;
}


int col_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c)
{
    if (C<=c) { fprintf(stderr,"error in col_inplace_c: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (iscolmajor)
    {
        size_t m = 2*c*R, n = 0;
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                if (m!=n) { cblas_ccopy((int)R,&X[m],1,&X[n],1); }
                m += 2*R*C; n += 2*R;
            }
        }
    }
    else if (c>0)
    {
        const size_t HS = H*S;
        if (HS==1) { cblas_ccopy((int)R,&X[2*c],(int)C,X,1); }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_ccopy((int)HS,&X[2*(c*HS+r*HS*C)],1,&X[2*r*HS],1);
            }
        }
    }

    return 0;
}


int col_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t c)
{
    if (C<=c) { fprintf(stderr,"error in col_inplace_z: C (ncols X) must be greater than c (col num to select)\n"); return 1; }

    if (iscolmajor)
    {
        size_t m = 2*c*R, n = 0;
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                if (m!=n) { cblas_zcopy((int)R,&X[m],1,&X[n],1); }
                m += 2*R*C; n += 2*R;
            }
        }
    }
    else if (c>0)
    {
        const size_t HS = H*S;
        if (HS==1) { cblas_zcopy((int)R,&X[2*c],(int)C,X,1); }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_zcopy((int)HS,&X[2*(c*HS+r*HS*C)],1,&X[2*r*HS],1);
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
