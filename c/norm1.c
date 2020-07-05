//Gets L1 norm (sum of absolute-values) for each row or col of X according to dim.
//For complex case, output is real.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int norm1_s (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int norm1_d (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int norm1_c (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int norm1_z (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);


int norm1_s (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    int r, c;

    if (R<1) { fprintf(stderr,"error in norm1_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in norm1_s: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                Y[c] = cblas_sasum((int)R,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                Y[c] = cblas_sasum((int)R,&X[c],(int)C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                Y[r] = cblas_sasum((int)C,&X[r],(int)R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                Y[r] = cblas_sasum((int)C,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in norm1_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int norm1_d (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    int r, c;

    if (R<1) { fprintf(stderr,"error in norm1_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in norm1_d: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                Y[c] = cblas_dasum((int)R,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                Y[c] = cblas_dasum((int)R,&X[c],(int)C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                Y[r] = cblas_dasum((int)C,&X[r],(int)R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                Y[r] = cblas_dasum((int)C,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in norm1_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int norm1_c (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    int r, c;

    if (R<1) { fprintf(stderr,"error in norm1_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in norm1_c: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                Y[c] = cblas_scasum((int)R,&X[2*c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                Y[c] = cblas_scasum((int)R,&X[2*c],(int)C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                Y[r] = cblas_scasum((int)C,&X[2*r],(int)R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                Y[r] = cblas_scasum((int)C,&X[2*r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in norm1_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int norm1_z (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    int r, c;

    if (R<1) { fprintf(stderr,"error in norm1_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in norm1_z: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                Y[c] = cblas_dzasum((int)R,&X[2*c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                Y[c] = cblas_dzasum((int)R,&X[2*c],(int)C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                Y[r] = cblas_dzasum((int)C,&X[2*r],(int)R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                Y[r] = cblas_dzasum((int)C,&X[2*r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in norm1_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
