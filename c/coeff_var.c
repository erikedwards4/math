//Gets the coefficient of variation (std/mean) of each row or col of X according to dim.
//For complex case, output is real.
//This works in place.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int coeff_var_s (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased);
int coeff_var_d (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased);


int coeff_var_s (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    const float o = 1.0f;
    float m, den;
    int r, c;

    if (R<1) { fprintf(stderr,"error in coeff_var_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in coeff_var_s: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = (float)R; } else { den = (float)(R-1); }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                m = cblas_sdot((int)R,&X[c*R],1,&o,0) / R;
                cblas_saxpy((int)R,-m,&o,0,&X[c*R],1);
                Y[c] = sqrtf(cblas_snrm2((int)R,&X[c*R],1)/den) / m;
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                m = cblas_sdot((int)R,&X[c],(int)C,&o,0) / R;
                cblas_saxpy((int)R,-m,&o,0,&X[c],(int)C);
                Y[c] = sqrtf(cblas_snrm2((int)R,&X[c],(int)C)/den) / m;
            }
        }
    }
    else if (dim==1)
    {
        if (biased) { den = (float)C; } else { den = (float)(C-1); }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                m = cblas_sdot((int)C,&X[r],(int)R,&o,0) / C;
                cblas_saxpy((int)C,-m,&o,0,&X[r],(int)R);
                Y[r] = sqrtf(cblas_snrm2((int)C,&X[r],(int)R)/den) / m;
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                m = cblas_sdot((int)C,&X[r*C],1,&o,0) / C;
                cblas_saxpy((int)C,-m,&o,0,&X[r*C],1);
                Y[r] = sqrtf(cblas_snrm2((int)C,&X[r*C],1)/den) / m;
            }
        }
    }
    else
    {
        fprintf(stderr,"error in coeff_var_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int coeff_var_d (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor, const char biased)
{
    const double o = 1.0;
    double m, den;
    int r, c;

    if (R<1) { fprintf(stderr,"error in coeff_var_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in coeff_var_d: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = (double)R; } else { den = (double)(R-1); }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                m = cblas_ddot((int)R,&X[c*R],1,&o,0) / R;
                cblas_daxpy((int)R,-m,&o,0,&X[c*R],1);
                Y[c] = sqrt(cblas_dnrm2((int)R,&X[c*R],1)/den) / m;
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                m = cblas_ddot((int)R,&X[c],(int)C,&o,0) / R;
                cblas_daxpy((int)R,-m,&o,0,&X[c],(int)C);
                Y[c] = sqrt(cblas_dnrm2((int)R,&X[c],(int)C)/den) / m;
            }
        }
    }
    else if (dim==1)
    {
        if (biased) { den = (double)C; } else { den = (double)(C-1); }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                m = cblas_ddot((int)C,&X[r],(int)R,&o,0) / C;
                cblas_daxpy((int)C,-m,&o,0,&X[r],(int)R);
                Y[r] = sqrt(cblas_dnrm2((int)C,&X[r],(int)R)/den) / m;
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                m = cblas_ddot((int)C,&X[r*C],1,&o,0) / C;
                cblas_daxpy((int)C,-m,&o,0,&X[r*C],1);
                Y[r] = sqrt(cblas_dnrm2((int)C,&X[r*C],1)/den) / m;
            }
        }
    }
    else
    {
        fprintf(stderr,"error in coeff_var_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
