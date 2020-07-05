//Joins 2 inputs X1, X2 into 1 output Y

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int join3_s (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t R2, const size_t R3, const size_t C1, const size_t C2, const size_t C3, const int dim, const char iscolmajor);
int join3_d (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t R2, const size_t R3, const size_t C1, const size_t C2, const size_t C3, const int dim, const char iscolmajor);
int join3_c (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t R2, const size_t R3, const size_t C1, const size_t C2, const size_t C3, const int dim, const char iscolmajor);
int join3_z (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t R2, const size_t R3, const size_t C1, const size_t C2, const size_t C3, const int dim, const char iscolmajor);


int join3_s (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t R2, const size_t R3, const size_t C1, const size_t C2, const size_t C3, const int dim, const char iscolmajor)
{
    const size_t R = R1+R2+R3, C = C1+C2+C3;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_scopy((int)R1,&X1[c*R1],1,&Y[c*R],1);
                cblas_scopy((int)R2,&X2[c*R2],1,&Y[c*R+R1],1);
                cblas_scopy((int)R3,&X3[c*R3],1,&Y[c*R+R1+R2],1);
            }
        }
        else
        {
            cblas_scopy((int)(R1*C1),X1,1,Y,1);
            cblas_scopy((int)(R2*C2),X2,1,&Y[R1*C1],1);
            cblas_scopy((int)(R3*C3),X3,1,&Y[R1*C1+R2*C2],1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_scopy((int)(R1*C1),X1,1,Y,1);
            cblas_scopy((int)(R2*C2),X2,1,&Y[R1*C1],1);
            cblas_scopy((int)(R3*C3),X3,1,&Y[R1*C1+R2*C2],1);
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_scopy((int)R1,&X1[r*C1],1,&Y[r*C],1);
                cblas_scopy((int)R2,&X2[r*C2],1,&Y[r*C+C1],1);
                cblas_scopy((int)R3,&X3[r*C3],1,&Y[r*C+C1+C2],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in join3_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int join3_d (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t R2, const size_t R3, const size_t C1, const size_t C2, const size_t C3, const int dim, const char iscolmajor)
{
    const size_t R = R1+R2+R3, C = C1+C2+C3;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_dcopy((int)R1,&X1[c*R1],1,&Y[c*R],1);
                cblas_dcopy((int)R2,&X2[c*R2],1,&Y[c*R+R1],1);
                cblas_dcopy((int)R3,&X3[c*R3],1,&Y[c*R+R1+R2],1);
            }
        }
        else
        {
            cblas_dcopy((int)(R1*C1),X1,1,Y,1);
            cblas_dcopy((int)(R2*C2),X2,1,&Y[R1*C1],1);
            cblas_dcopy((int)(R3*C3),X3,1,&Y[R1*C1+R2*C2],1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_dcopy((int)(R1*C1),X1,1,Y,1);
            cblas_dcopy((int)(R2*C2),X2,1,&Y[R1*C1],1);
            cblas_dcopy((int)(R3*C3),X3,1,&Y[R1*C1+R2*C2],1);
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_dcopy((int)R1,&X1[r*C1],1,&Y[r*C],1);
                cblas_dcopy((int)R2,&X2[r*C2],1,&Y[r*C+C1],1);
                cblas_dcopy((int)R3,&X3[r*C3],1,&Y[r*C+C1+C2],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in join3_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int join3_c (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t R2, const size_t R3, const size_t C1, const size_t C2, const size_t C3, const int dim, const char iscolmajor)
{
    const size_t R = R1+R2+R3, C = C1+C2+C3;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_ccopy((int)R1,&X1[2*c*R1],1,&Y[2*c*R],1);
                cblas_ccopy((int)R2,&X2[2*c*R2],1,&Y[2*(c*R+R1)],1);
                cblas_ccopy((int)R3,&X3[2*c*R3],1,&Y[2*(c*R+R1+R2)],1);
            }
        }
        else
        {
            cblas_ccopy((int)(R1*C1),X1,1,Y,1);
            cblas_ccopy((int)(R2*C2),X2,1,&Y[2*R1*C1],1);
            cblas_ccopy((int)(R3*C3),X3,1,&Y[2*(R1*C1+R2*C2)],1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_ccopy((int)(R1*C1),X1,1,Y,1);
            cblas_ccopy((int)(R2*C2),X2,1,&Y[2*R1*C1],1);
            cblas_ccopy((int)(R3*C3),X3,1,&Y[2*(R1*C1+R2*C2)],1);
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_ccopy((int)R1,&X1[2*r*C1],1,&Y[2*r*C],1);
                cblas_ccopy((int)R2,&X2[2*r*C2],1,&Y[2*(r*C+C1)],1);
                cblas_ccopy((int)R3,&X3[2*r*C3],1,&Y[2*(r*C+C1+C2)],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in join3_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int join3_z (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t R2, const size_t R3, const size_t C1, const size_t C2, const size_t C3, const int dim, const char iscolmajor)
{
    const size_t R = R1+R2+R3, C = C1+C2+C3;

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_zcopy((int)R1,&X1[2*c*R1],1,&Y[2*c*R],1);
                cblas_zcopy((int)R2,&X2[2*c*R2],1,&Y[2*(c*R+R1)],1);
                cblas_zcopy((int)R3,&X3[2*c*R3],1,&Y[2*(c*R+R1+R2)],1);
            }
        }
        else
        {
            cblas_zcopy((int)(R1*C1),X1,1,Y,1);
            cblas_zcopy((int)(R2*C2),X2,1,&Y[2*R1*C1],1);
            cblas_zcopy((int)(R3*C3),X3,1,&Y[2*(R1*C1+R2*C2)],1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_zcopy((int)(R1*C1),X1,1,Y,1);
            cblas_zcopy((int)(R2*C2),X2,1,&Y[2*R1*C1],1);
            cblas_zcopy((int)(R3*C3),X3,1,&Y[2*(R1*C1+R2*C2)],1);
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_zcopy((int)R1,&X1[2*r*C1],1,&Y[2*r*C],1);
                cblas_zcopy((int)R2,&X2[2*r*C2],1,&Y[2*(r*C+C1)],1);
                cblas_zcopy((int)R3,&X3[2*r*C3],1,&Y[2*(r*C+C1+C2)],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in join3_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
