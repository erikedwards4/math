//Gets one row of input X

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int row_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r);
int row_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r);
int row_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r);
int row_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r);

int row_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r);
int row_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r);
int row_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r);
int row_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r);


int row_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_s: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        X += r;
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, X+=R*C, Y+=C)
            {
                cblas_scopy((int)C,X,(int)R,Y,1);
            }
        }
    }
    else if (r>0)
    {
        cblas_scopy((int)(H*S*C),&X[r*H*S*C],1,Y,1);
    }

    return 0;
}


int row_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_d: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        X += r;
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, X+=R*C, Y+=C)
            {
                cblas_dcopy((int)C,X,(int)R,Y,1);
            }
        }
    }
    else if (r>0)
    {
        cblas_dcopy((int)(H*S*C),&X[r*H*S*C],1,Y,1);
    }

    return 0;
}


int row_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_c: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        X += 2*r;
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, X+=2*R*C, Y+=2*C)
            {
                cblas_ccopy((int)C,X,(int)R,Y,1);
            }
        }
    }
    else if (r>0)
    {
        cblas_ccopy((int)(H*S*C),&X[2*r*H*S*C],1,Y,1);
    }

    return 0;
}


int row_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_z: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        X += 2*r;
        for (size_t h=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++, X+=2*R*C, Y+=2*C)
            {
                cblas_zcopy((int)C,X,(int)R,Y,1);
            }
        }
    }
    else if (r>0)
    {
        cblas_zcopy((int)(H*S*C),&X[2*r*H*S*C],1,Y,1);
    }

    return 0;
}


int row_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_inplace_s: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        for (size_t h=0, n=0, m=r; h<H; h++)
        {
            for (size_t s=0; s<S; s++, m+=R*C, n+=C)
            {
                cblas_scopy((int)C,&X[m],(int)R,&X[n],1);
            }
        }
    }
    else if (r>0)
    {
        cblas_scopy((int)(H*S*C),&X[r*H*S*C],1,X,1);
    }

    return 0;
}


int row_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_inplace_d: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        for (size_t h=0, n=0, m=r; h<H; h++)
        {
            for (size_t s=0; s<S; s++, m+=R*C, n+=C)
            {
                cblas_dcopy((int)C,&X[m],(int)R,&X[n],1);
            }
        }
    }
    else if (r>0)
    {
        cblas_dcopy((int)(H*S*C),&X[r*H*S*C],1,X,1);
    }

    return 0;
}


int row_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_inplace_c: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        for (size_t h=0, n=0, m=2*r; h<H; h++)
        {
            for (size_t s=0; s<S; s++, m+=2*R*C, n+=2*C)
            {
                cblas_ccopy((int)C,&X[m],(int)R,&X[n],1);
            }
        }
    }
    else if (r>0)
    {
        cblas_ccopy((int)(H*S*C),&X[2*r*H*S*C],1,X,1);
    }

    return 0;
}


int row_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t r)
{
    if (r>R) { fprintf(stderr,"error in row_inplace_z: R (nrows X) must be greater than r (row num to select)\n"); return 1; }

    if (iscolmajor)
    {
        for (size_t h=0, n=0, m=2*r; h<H; h++)
        {
            for (size_t s=0; s<S; s++, m+=2*R*C, n+=2*C)
            {
                cblas_zcopy((int)C,&X[m],(int)R,&X[n],1);
            }
        }
    }
    else if (r>0)
    {
        cblas_zcopy((int)(H*S*C),&X[2*r*H*S*C],1,X,1);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
