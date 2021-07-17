//Linear algebra function.
//Multiplies a square matrix by itself: Y = X*X.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int matmul1_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor);
int matmul1_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor);
int matmul1_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor);
int matmul1_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor);


int matmul1_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R!=C) { fprintf(stderr,"error in matmul1_s: X must be a square matrix \n"); return 1; }
    
    const size_t N = R*C;
    float sm2;

    if (N==0u) {}
    else if (N<1100u)
    {
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                for (size_t r1=0u; r1<R; ++r1, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t r2=0u; r2<R; ++r2)
                    {
                        //sm2 += *(X+r2*R+r1) * *(X+c*R+r2);
                        sm2 = fmaf(*(X+r2*R+r1),*(X+c*R+r2),sm2);
                    }
                    *Y = sm2;
                }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t c2=0u; c2<C; ++c2, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t c1=0u; c1<C; ++c1)
                    {
                        //sm2 += *(X+r*C+c1) * *(X+c1*R+c2);
                        sm2 = fmaf(*(X+r*C+c1),*(X+c1*R+c2),sm2);
                    }
                    *Y = sm2;
                }
            }
        }
    }
    else
    {
        if (iscolmajor)
        {
            cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,(int)R,1.0f,X,(int)R,X,(int)R,0.0f,Y,(int)R);
        }
        else
        {
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,(int)C,1.0f,X,(int)C,X,(int)C,0.0f,Y,(int)C);
        }
    }

    return 0;
}


int matmul1_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R!=C) { fprintf(stderr,"error in matmul1_d: X must be a square matrix \n"); return 1; }
    
    const size_t N = R*C;
    double sm2;

    if (N==0u) {}
    else if (N<1100u)
    {
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                for (size_t r1=0u; r1<R; ++r1, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t r2=0u; r2<R; ++r2)
                    {
                        sm2 = fma(*(X+r2*R+r1),*(X+c*R+r2),sm2);
                    }
                    *Y = sm2;
                }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t c2=0u; c2<C; ++c2, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t c1=0u; c1<C; ++c1)
                    {
                        sm2 = fma(*(X+r*C+c1),*(X+c1*R+c2),sm2);
                    }
                    *Y = sm2;
                }
            }
        }
    }
    else
    {
        if (iscolmajor)
        {
            cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,(int)R,1.0,X,(int)R,X,(int)R,0.0,Y,(int)R);
        }
        else
        {
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,(int)C,1.0,X,(int)C,X,(int)C,0.0,Y,(int)C);
        }
    }

    return 0;
}


int matmul1_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R!=C) { fprintf(stderr,"error in matmul1_c: X must be a square matrix \n"); return 1; }
    
    const size_t N = R*C;
    float x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0u) {}
    else if (N<1100u)
    {
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                for (size_t r1=0u; r1<R; ++r1, ++Y)
                {
                    sm2r = sm2i = 0.0f;
                    for (size_t r2=0u; r2<R; ++r2)
                    {
                        x1r = *(X+2u*(r2*R+r1)); x1i = *(X+2u*(r2*R+r1)+1u);
                        x2r = *(X+2u*(c*R+r2)); x2i = *(X+2u*(c*R+r2)+1u);
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t c2=0u; c2<C; ++c2, ++Y)
                {
                    sm2r = sm2i = 0.0f;
                    for (size_t c1=0u; c1<C; ++c1)
                    {
                        x1r = *(X+2u*(r*C+c1)); x1i = *(X+2u*(r*C+c1)+1u);
                        x2r = *(X+2u*(c1*R+c2)); x2i = *(X+2u*(c1*R+c2)+1u);
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                }
            }
        }
    }
    else
    {
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
        if (iscolmajor)
        {
            cblas_cgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,(int)R,o,X,(int)R,X,(int)R,z,Y,(int)R);
        }
        else
        {
            cblas_cgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,(int)C,o,X,(int)C,X,(int)C,z,Y,(int)C);
        }
    }

    return 0;
}


int matmul1_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R!=C) { fprintf(stderr,"error in matmul1_z: X must be a square matrix \n"); return 1; }
    
    const size_t N = R*C;
    double x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0u) {}
    else if (N<1100u)
    {
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                for (size_t r1=0u; r1<R; ++r1, ++Y)
                {
                    sm2r = sm2i = 0.0;
                    for (size_t r2=0u; r2<R; ++r2)
                    {
                        x1r = *(X+2u*(r2*R+r1)); x1i = *(X+2u*(r2*R+r1)+1u);
                        x2r = *(X+2u*(c*R+r2)); x2i = *(X+2u*(c*R+r2)+1u);
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                for (size_t c2=0u; c2<C; ++c2, ++Y)
                {
                    sm2r = sm2i = 0.0;
                    for (size_t c1=0u; c1<C; ++c1)
                    {
                        x1r = *(X+2u*(r*C+c1)); x1i = *(X+2u*(r*C+c1)+1u);
                        x2r = *(X+2u*(c1*R+c2)); x2i = *(X+2u*(c1*R+c2)+1u);
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                }
            }
        }
    }
    else
    {
        const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
        if (iscolmajor)
        {
            cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,(int)R,o,X,(int)R,X,(int)R,z,Y,(int)R);
        }
        else
        {
            cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R,(int)C,(int)C,o,X,(int)C,X,(int)C,z,Y,(int)C);
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
