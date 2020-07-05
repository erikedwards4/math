//Splits 1 input X into 3 outputs Y1, Y2, Y3

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int split3_s (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const int dim, const char iscolmajor);
int split3_d (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const int dim, const char iscolmajor);
int split3_c (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const int dim, const char iscolmajor);
int split3_z (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const int dim, const char iscolmajor);


int split3_s (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const int dim, const char iscolmajor)
{
    const size_t R3 = R/3, C3 = C/3;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_scopy((int)R3,&X[c*R],1,&Y1[c*R3],1);
                cblas_scopy((int)R3,&X[c*R+R3],1,&Y2[c*R3],1);
                cblas_scopy((int)R3,&X[c*R+2*R3],1,&Y3[c*R3],1);
            }
        }
        else
        {
            cblas_scopy((int)(R*C3),X,1,Y1,1);
            cblas_scopy((int)(R*C3),&X[R*C3],1,Y2,1);
            cblas_scopy((int)(R*C3),&X[2*R*C3],1,Y3,1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_scopy((int)(R*C3),X,1,Y1,1);
            cblas_scopy((int)(R*C3),&X[R*C3],1,Y2,1);
            cblas_scopy((int)(R*C3),&X[2*R*C3],1,Y3,1);
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_scopy((int)C3,&X[r*C],1,&Y1[r*C3],1);
                cblas_scopy((int)C3,&X[r*C+C3],1,&Y2[r*C3],1);
                cblas_scopy((int)C3,&X[r*C+2*C3],1,&Y3[r*C3],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in split3_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int split3_d (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const int dim, const char iscolmajor)
{
    const size_t R3 = R/3, C3 = C/3;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_dcopy((int)R3,&X[c*R],1,&Y1[c*R3],1);
                cblas_dcopy((int)R3,&X[c*R+R3],1,&Y2[c*R3],1);
                cblas_dcopy((int)R3,&X[c*R+2*R3],1,&Y3[c*R3],1);
            }
        }
        else
        {
            cblas_dcopy((int)(R*C3),X,1,Y1,1);
            cblas_dcopy((int)(R*C3),&X[R*C3],1,Y2,1);
            cblas_dcopy((int)(R*C3),&X[2*R*C3],1,Y3,1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_dcopy((int)(R*C3),X,1,Y1,1);
            cblas_dcopy((int)(R*C3),&X[R*C3],1,Y2,1);
            cblas_dcopy((int)(R*C3),&X[2*R*C3],1,Y3,1);
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_dcopy((int)C3,&X[r*C],1,&Y1[r*C3],1);
                cblas_dcopy((int)C3,&X[r*C+C3],1,&Y2[r*C3],1);
                cblas_dcopy((int)C3,&X[r*C+2*C3],1,&Y3[r*C3],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in split3_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int split3_c (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const int dim, const char iscolmajor)
{
    const size_t R3 = R/3, C3 = C/3;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_ccopy((int)R3,&X[2*c*R],1,&Y1[2*c*R3],1);
                cblas_ccopy((int)R3,&X[2*(c*R+R3)],1,&Y2[2*c*R3],1);
                cblas_ccopy((int)R3,&X[2*(c*R+2*R3)],1,&Y3[2*c*R3],1);
            }
        }
        else
        {
            cblas_ccopy((int)(R*C3),X,1,Y1,1);
            cblas_ccopy((int)(R*C3),&X[2*R*C3],1,Y2,1);
            cblas_ccopy((int)(R*C3),&X[4*R*C3],1,Y3,1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_ccopy((int)(R*C3),X,1,Y1,1);
            cblas_ccopy((int)(R*C3),&X[2*R*C3],1,Y2,1);
            cblas_ccopy((int)(R*C3),&X[4*R*C3],1,Y3,1);
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_ccopy((int)C3,&X[2*r*C],1,&Y1[2*r*C3],1);
                cblas_ccopy((int)C3,&X[2*(r*C+C3)],1,&Y2[2*r*C3],1);
                cblas_ccopy((int)C3,&X[2*(r*C+2*C3)],1,&Y3[2*r*C3],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in split3_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int split3_z (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const int dim, const char iscolmajor)
{
    const size_t R3 = R/3, C3 = C/3;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_zcopy((int)R3,&X[2*c*R],1,&Y1[2*c*R3],1);
                cblas_zcopy((int)R3,&X[2*(c*R+R3)],1,&Y2[2*c*R3],1);
                cblas_zcopy((int)R3,&X[2*(c*R+2*R3)],1,&Y3[2*c*R3],1);
            }
        }
        else
        {
            cblas_zcopy((int)(R*C3),X,1,Y1,1);
            cblas_zcopy((int)(R*C3),&X[2*R*C3],1,Y2,1);
            cblas_zcopy((int)(R*C3),&X[4*R*C3],1,Y3,1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_zcopy((int)(R*C3),X,1,Y1,1);
            cblas_zcopy((int)(R*C3),&X[2*R*C3],1,Y2,1);
            cblas_zcopy((int)(R*C3),&X[4*R*C3],1,Y3,1);
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_zcopy((int)C3,&X[2*r*C],1,&Y1[2*r*C3],1);
                cblas_zcopy((int)C3,&X[2*(r*C+C3)],1,&Y2[2*r*C3],1);
                cblas_zcopy((int)C3,&X[2*(r*C+2*C3)],1,&Y3[2*r*C3],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in split3_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
