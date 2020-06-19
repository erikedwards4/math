//Gets skewness of each row or col of X according to dim.
//For complex case, output is real.
//This works in place.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int skewness_s (float *Y, float *X, const int R, const int C, const int dim, const char iscolmajor, const char biased);
int skewness_d (double *Y, double *X, const int R, const int C, const int dim, const char iscolmajor, const char biased);
int skewness_c (float *Y, float *X, const int R, const int C, const int dim, const char iscolmajor, const char biased);
int skewness_z (double *Y, double *X, const int R, const int C, const int dim, const char iscolmajor, const char biased);


int skewness_s (float *Y, float *X, const int R, const int C, const int dim, const char iscolmajor, const char biased)
{
    const float o = 1.0f;
    float m, s2, s3;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in skewness_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in skewness_s: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = (float)R; } else { den = (float)(R-1); }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                m = -cblas_sdot(R,&X[c*R],1,&o,0) / R;
                cblas_saxpy(R,m,&o,0,&X[c*R],1);
                s2 = s3 = 0.0f;
                for (n=c*R; n<(c+1)*R; n++) { x2 = X[n]*X[n]; s2 += x2; x2 *= X[n]; s3 += x2; }
                Y[c] = adj * s3 / (s2*sqrtf(s2));
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                m = -cblas_sdot(R,&X[c],C,&o,0) / R;
                cblas_saxpy(R,m,&o,0,&X[c],C);
                Y[c] = cblas_snrm2(R,&X[c],C) / den;
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
                m = -cblas_sdot(C,&X[r],R,&o,0) / C;
                cblas_saxpy(C,m,&o,0,&X[r],R);
                Y[r] = cblas_snrm2(C,&X[r],R) / den;
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                m = -cblas_sdot(C,&X[r*C],1,&o,0) / C;
                cblas_saxpy(C,m,&o,0,&X[r*C],1);
                Y[r] = cblas_snrm2(C,&X[r*C],1) / den;
            }
        }
    }
    else
    {
        fprintf(stderr,"error in skewness_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}
