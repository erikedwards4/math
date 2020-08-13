//Linear algebra function.
//Multiplies a square matrix by itself transposed: Y = X*X'.

//For complex case, this uses the Hermitian (conjugate) transpose.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int matmul1t_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor);
int matmul1t_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor);
int matmul1t_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor);
int matmul1t_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor);


int matmul1t_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor)
{
    const size_t N = R*C;
    float sm2;

    if (N==0) {}
    else if (N<1100)
    {
        if (iscolmajor)
        {
            for (size_t r2=0; r2<R; ++r2)
            {
                for (size_t r1=0; r1<R; ++r1, X-=N, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t c=0; c<C; ++c, X+=R)
                    {
                        sm2 = fmaf(*(X+r1),*(X+r2),sm2);
                    }
                    *Y = sm2;
                }
            }
        }
        else
        {
            for (size_t r1=0; r1<R; ++r1)
            {
                for (size_t r2=0; r2<R; ++r2, X-=C, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t c=0; c<C; ++c, ++X)
                    {
                        sm2 = fmaf(*(X+r1*C),*(X+r2*C),sm2);
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
            cblas_sgemm(CblasColMajor,CblasNoTrans,CblasTrans,(int)R,(int)C,(int)R,1.0f,X,(int)R,X,(int)R,0.0f,Y,(int)R);
        }
        else
        {
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)R,(int)C,(int)C,1.0f,X,(int)C,X,(int)C,0.0f,Y,(int)C);
        }
    }

    return 0;
}


int matmul1t_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor)
{
    const size_t N = R*C;
    double sm2;

    if (N==0) {}
    else if (N<1100)
    {
        if (iscolmajor)
        {
            for (size_t r2=0; r2<R; ++r2)
            {
                for (size_t r1=0; r1<R; ++r1, X-=N, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t c=0; c<C; ++c, X+=R)
                    {
                        sm2 = fma(*(X+r1),*(X+r2),sm2);
                    }
                    *Y = sm2;
                }
            }
        }
        else
        {
            for (size_t r1=0; r1<R; ++r1)
            {
                for (size_t r2=0; r2<R; ++r2, X-=C, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t c=0; c<C; ++c, ++X)
                    {
                        sm2 = fma(*(X+r1*C),*(X+r2*C),sm2);
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
            cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,(int)R,(int)C,(int)R,1.0,X,(int)R,X,(int)R,1.0,Y,(int)R);
        }
        else
        {
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)R,(int)C,(int)C,1.0,X,(int)C,X,(int)C,1.0,Y,(int)C);
        }
    }

    return 0;
}


int matmul1t_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor)
{
    const size_t N = R*C;
    float x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0) {}
    else if (N<1100)
    {
        if (iscolmajor)
        {
            for (size_t r2=0; r2<R; ++r2)
            {
                for (size_t r1=0; r1<R; ++r1, X-=2*N, ++Y)
                {
                    sm2r = sm2i = 0.0f;
                    for (size_t c=0; c<C; ++c, X+=2*R)
                    {
                        x1r = *(X+2*r1); x1i = *(X+2*r1+1);
                        x2r = *(X+2*r2); x2i = -*(X+2*r2+1);
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                }
            }
        }
        else
        {
            for (size_t r1=0; r1<R; ++r1)
            {
                for (size_t r2=0; r2<R; ++r2, X-=2*C, ++Y)
                {
                    sm2r = sm2i = 0.0f;
                    for (size_t c=0; c<C; ++c, X+=2)
                    {
                        x1r = *(X+2*r1*C); x1i = *(X+2*r1*C+1);
                        x2r = *(X+2*r2*C); x2i = -*(X+2*r2*C+1);
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
        const float o[2] = {1.0f,0.0f};
        if (iscolmajor)
        {
            cblas_cgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,(int)R,(int)C,(int)R,o,X,(int)R,X,(int)R,o,Y,(int)R);
        }
        else
        {
            cblas_cgemm(CblasRowMajor,CblasNoTrans,CblasConjTrans,(int)R,(int)C,(int)C,o,X,(int)C,X,(int)C,o,Y,(int)C);
        }
    }

    return 0;
}


int matmul1t_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor)
{
    const size_t N = R*C;
    double x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0) {}
    else if (N<1100)
    {
        if (iscolmajor)
        {
            for (size_t r2=0; r2<R; ++r2)
            {
                for (size_t r1=0; r1<R; ++r1, X-=2*N, ++Y)
                {
                    sm2r = sm2i = 0.0;
                    for (size_t c=0; c<C; ++c, X+=2*R)
                    {
                        x1r = *(X+2*r1); x1i = *(X+2*r1+1);
                        x2r = *(X+2*r2); x2i = -*(X+2*r2+1);
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                }
            }
        }
        else
        {
            for (size_t r1=0; r1<R; ++r1)
            {
                for (size_t r2=0; r2<R; ++r2, X-=2*C, ++Y)
                {
                    sm2r = sm2i = 0.0;
                    for (size_t c=0; c<C; ++c, X+=2)
                    {
                        x1r = *(X+2*r1*C); x1i = *(X+2*r1*C+1);
                        x2r = *(X+2*r2*C); x2i = -*(X+2*r2*C+1);
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
        const double o[2] = {1.0,0.0};
        if (iscolmajor)
        {
            cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,(int)R,(int)C,(int)R,o,X,(int)R,X,(int)R,o,Y,(int)R);
        }
        else
        {
            cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasConjTrans,(int)R,(int)C,(int)C,o,X,(int)C,X,(int)C,o,Y,(int)C);
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
