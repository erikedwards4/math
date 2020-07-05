//Splits 1 input X into 2 outputs Y1, Y2

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int split2_s (float *Y1, float *Y2, const float *X, const size_t R, const size_t C, const int dim, const char iscolmajor);
int split2_d (double *Y1, double *Y2, const double *X, const size_t R, const size_t C, const int dim, const char iscolmajor);
int split2_c (float *Y1, float *Y2, const float *X, const size_t R, const size_t C, const int dim, const char iscolmajor);
int split2_z (double *Y1, double *Y2, const double *X, const size_t R, const size_t C, const int dim, const char iscolmajor);


int split2_s (float *Y1, float *Y2, const float *X, const size_t R, const size_t C, const int dim, const char iscolmajor)
{
    const size_t R2 = R/2, C2 = C/2;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_scopy((int)R2,&X[c*R],1,&Y1[c*R2],1);
                cblas_scopy((int)R2,&X[c*R+R2],1,&Y2[c*R2],1);
            }
        }
        else
        {
            cblas_scopy((int)(R*C2),X,1,Y1,1);
            cblas_scopy((int)(R*C2),&X[R*C2],1,Y2,1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_scopy((int)(R*C2),X,1,Y1,1);
            cblas_scopy((int)(R*C2),&X[R*C2],1,Y2,1);
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_scopy((int)C2,&X[r*C],1,&Y1[r*C2],1);
                cblas_scopy((int)C2,&X[r*C+C2],1,&Y2[r*C2],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in split2_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int split2_d (double *Y1, double *Y2, const double *X, const size_t R, const size_t C, const int dim, const char iscolmajor)
{
    const size_t R2 = R/2, C2 = C/2;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_dcopy((int)R2,&X[c*R],1,&Y1[c*R2],1);
                cblas_dcopy((int)R2,&X[c*R+R2],1,&Y2[c*R2],1);
            }
        }
        else
        {
            cblas_dcopy((int)(R*C2),X,1,Y1,1);
            cblas_dcopy((int)(R*C2),&X[R*C2],1,Y2,1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_dcopy((int)(R*C2),X,1,Y1,1);
            cblas_dcopy((int)(R*C2),&X[R*C2],1,Y2,1);
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_dcopy((int)C2,&X[r*C],1,&Y1[r*C2],1);
                cblas_dcopy((int)C2,&X[r*C+C2],1,&Y2[r*C2],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in split2_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int split2_c (float *Y1, float *Y2, const float *X, const size_t R, const size_t C, const int dim, const char iscolmajor)
{
    const size_t R2 = R/2, C2 = C/2;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_ccopy((int)R2,&X[2*c*R],1,&Y1[c*R],1);
                cblas_ccopy((int)R2,&X[2*c*R+R],1,&Y2[c*R],1);
            }
        }
        else
        {
            cblas_ccopy((int)(R*C2),X,1,Y1,1);
            cblas_ccopy((int)(R*C2),&X[R*C],1,Y2,1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_ccopy((int)(R*C2),X,1,Y1,1);
            cblas_ccopy((int)(R*C2),&X[R*C],1,Y2,1);
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_ccopy((int)C2,&X[2*r*C],1,&Y1[r*C],1);
                cblas_ccopy((int)C2,&X[2*r*C+C],1,&Y2[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in split2_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int split2_z (double *Y1, double *Y2, const double *X, const size_t R, const size_t C, const int dim, const char iscolmajor)
{
    const size_t R2 = R/2, C2 = C/2;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_zcopy((int)R2,&X[2*c*R],1,&Y1[c*R],1);
                cblas_zcopy((int)R2,&X[2*c*R+R],1,&Y2[c*R],1);
            }
        }
        else
        {
            cblas_zcopy((int)(R*C2),X,1,Y1,1);
            cblas_zcopy((int)(R*C2),&X[R*C],1,Y2,1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_zcopy((int)(R*C2),X,1,Y1,1);
            cblas_zcopy((int)(R*C2),&X[R*C],1,Y2,1);
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_zcopy((int)C2,&X[2*r*C],1,&Y1[r*C],1);
                cblas_zcopy((int)C2,&X[2*r*C+C],1,&Y2[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in split2_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
