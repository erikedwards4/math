//Linear algebra function.
//Matrix multiplication of 2 input matrices X1 and X2,
//with the 2nd matrix transposed: Y = X1*X2'.

//For complex case, this uses the Hermitian (conjugate) transpose.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int matmul2t_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const char iscolmajor, const char tr);
int matmul2t_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const char iscolmajor, const char tr);
int matmul2t_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const char iscolmajor, const char tr);
int matmul2t_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const char iscolmajor, const char tr);


int matmul2t_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const char iscolmajor, const char tr)
{
    const size_t N = (tr) ? C1*C2 : R1*R2;
    float sm2;

    if (N==0) {}
    else if (N<1100)
    {
        const size_t N1 = R1*C1, N2 = R2*C2;
        if (tr)
        {
            if (iscolmajor)
            {
                for (size_t c2=0; c2<C2; ++c2)
                {
                    for (size_t c1=0; c1<C1; ++c1, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t r2=0; r2<R2; ++r2, ++X1, ++X2) { sm2 = fmaf(*X1,*X2,sm2); }
                        *Y = sm2;
                        if (c1<C1-1) { X2 -= R2; }
                    }
                    if (c2<C2-1) { X1 -= N1; }
                }
            }
            else
            {
                for (size_t c1=0; c1<C1; ++c1)
                {
                    for (size_t c2=0; c2<C2; ++c2, X1-=N1, X2-=N2-1, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t r1=0; r1<R1; ++r1, X1+=C1, X2+=C2) { sm2 = fmaf(*X1,*X2,sm2); }
                        *Y = sm2;
                    }
                    if (c1<C1-1) { ++X1; X2 -= C2; }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (size_t r2=0; r2<R2; ++r2)
                {
                    for (size_t r1=0; r1<R1; ++r1, X1-=N1-1, X2-=N2, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t c2=0; c2<C2; ++c2, X1+=R1, X2+=R2) { sm2 = fmaf(*X1,*X2,sm2); }
                        *Y = sm2;
                    }
                    if (r2<R2-1) { X1 -= R1; ++X2; }
                }
            }
            else
            {
                for (size_t r1=0; r1<R1; ++r1)
                {
                    for (size_t r2=0; r2<R2; ++r2, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t c1=0; c1<C1; ++c1, ++X1, ++X2) { sm2 = fmaf(*X1,*X2,sm2); }
                        *Y = sm2;
                        if (r2<R2-1) { X1 -= C1; }
                    }
                    if (r1<R1-1) { X2 -= N2; }
                }
            }
        }
    }
    else
    {
        if (tr)
        {
            if (iscolmajor)
            {
                cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,(int)C1,(int)C2,(int)R1,1.0f,X1,(int)R1,X2,(int)R1,0.0f,Y,(int)C1);
            }
            else
            {
                cblas_sgemm(CblasRowMajor,CblasTrans,CblasNoTrans,(int)C1,(int)C2,(int)R1,1.0f,X1,(int)C1,X2,(int)C2,0.0f,Y,(int)C2);
            }
        }
        else
        {
            if (iscolmajor)
            {
                cblas_sgemm(CblasColMajor,CblasNoTrans,CblasTrans,(int)R1,(int)R2,(int)C1,1.0f,X1,(int)R1,X2,(int)R2,0.0f,Y,(int)R1);
            }
            else
            {
                cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)R1,(int)R2,(int)C1,1.0f,X1,(int)C2,X2,(int)C2,0.0f,Y,(int)R2);
            }
        }
    }

    return 0;
}


int matmul2t_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const char iscolmajor, const char tr)
{
    const size_t N = (tr) ? C1*C2 : R1*R2;
    double sm2;

    if (N==0) {}
    else if (N<1100)
    {
        const size_t N1 = R1*C1, N2 = R2*C2;
        if (tr)
        {
            if (iscolmajor)
            {
                for (size_t c2=0; c2<C2; ++c2)
                {
                    for (size_t c1=0; c1<C1; ++c1, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t r2=0; r2<R2; ++r2, ++X1, ++X2) { sm2 = fma(*X1,*X2,sm2); }
                        *Y = sm2;
                        if (c1<C1-1) { X2 -= R2; }
                    }
                    if (c2<C2-1) { X1 -= N1; }
                }
            }
            else
            {
                for (size_t c1=0; c1<C1; ++c1)
                {
                    for (size_t c2=0; c2<C2; ++c2, X1-=N1, X2-=N2-1, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t r1=0; r1<R1; ++r1, X1+=C1, X2+=C2) { sm2 = fma(*X1,*X2,sm2); }
                        *Y = sm2;
                    }
                    if (c1<C1-1) { ++X1; X2 -= C2; }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (size_t r2=0; r2<R2; ++r2)
                {
                    for (size_t r1=0; r1<R1; ++r1, X1-=N1-1, X2-=N2, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t c2=0; c2<C2; ++c2, X1+=R1, X2+=R2) { sm2 = fma(*X1,*X2,sm2); }
                        *Y = sm2;
                    }
                    if (r2<R2-1) { X1 -= R1; ++X2; }
                }
            }
            else
            {
                for (size_t r1=0; r1<R1; ++r1)
                {
                    for (size_t r2=0; r2<R2; ++r2, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t c1=0; c1<C1; ++c1, ++X1, ++X2) { sm2 = fma(*X1,*X2,sm2); }
                        *Y = sm2;
                        if (r2<R2-1) { X1 -= C1; }
                    }
                    if (r1<R1-1) { X2 -= N2; }
                }
            }
        }
    }
    else
    {
        if (tr)
        {
            if (iscolmajor)
            {
                cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,(int)C1,(int)C2,(int)R1,1.0,X1,(int)R1,X2,(int)R1,0.0,Y,(int)C1);
            }
            else
            {
                cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,(int)C1,(int)C2,(int)R1,1.0,X1,(int)C1,X2,(int)C2,0.0,Y,(int)C2);
            }
        }
        else
        {
            if (iscolmajor)
            {
                cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,(int)R1,(int)R2,(int)C1,1.0,X1,(int)R1,X2,(int)R2,0.0,Y,(int)R1);
            }
            else
            {
                cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,(int)R1,(int)R2,(int)C1,1.0,X1,(int)C2,X2,(int)C2,0.0,Y,(int)R2);
            }
        }
    }

    return 0;
}


int matmul2t_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const char iscolmajor, const char tr)
{
    const size_t N = (tr) ? C1*C2 : R1*R2;
    float x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0) {}
    else if (N<1100)
    {
        const size_t N1 = R1*C1, N2 = R2*C2;
        if (tr)
        {
            if (iscolmajor)
            {
                for (size_t c2=0; c2<C2; ++c2)
                {
                    for (size_t c1=0; c1<C1; ++c1, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t r2=0; r2<R2; ++r2, ++X1, ++X2)
                        {
                            x1r = *X1; x1i = -*++X1;
                            x2r = *X2; x2i = *++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (c1<C1-1) { X2 -= 2*R2; }
                    }
                    if (c2<C2-1) { X1 -= 2*N1; }
                }
            }
            else
            {
                for (size_t c1=0; c1<C1; ++c1)
                {
                    for (size_t c2=0; c2<C2; ++c2, X1-=2*N1, X2-=2*N2-2, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t r1=0; r1<R1; ++r1, X1+=2*C1-1, X2+=2*C2-1)
                        {
                            x1r = *X1; x1i = -*++X1;
                            x2r = *X2; x2i = *++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                    }
                    if (c1<C1-1) { X1 += 2; X2 -= 2*C2; }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (size_t r2=0; r2<R2; ++r2)
                {
                    for (size_t r1=0; r1<R1; ++r1, X1-=2*N1-2, X2-=2*N2, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t c2=0; c2<C2; ++c2, X1+=2*R1-1, X2+=2*R2-1)
                        {
                            x1r = *X1; x1i = *++X1;
                            x2r = *X2; x2i = -*++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                    }
                    if (r2<R2-1) { X1 -= 2*R1; X2 += 2; }
                }
            }
            else
            {
                for (size_t r1=0; r1<R1; ++r1)
                {
                    for (size_t r2=0; r2<R2; ++r2, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t c1=0; c1<C1; ++c1, ++X1, ++X2)
                        {
                            x1r = *X1; x1i = *++X1;
                            x2r = *X2; x2i = -*++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (r2<R2-1) { X1 -= 2*C1; }
                    }
                    if (r1<R1-1) { X2 -= 2*N2; }
                }
            }
        }
    }
    else
    {
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
        if (tr)
        {
            if (iscolmajor)
            {
                cblas_cgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,(int)C1,(int)C2,(int)R1,o,X1,(int)R1,X2,(int)R1,z,Y,(int)C1);
            }
            else
            {
                cblas_cgemm(CblasRowMajor,CblasConjTrans,CblasNoTrans,(int)C1,(int)C2,(int)R1,o,X1,(int)C1,X2,(int)C2,z,Y,(int)C2);
            }
        }
        else
        {
            if (iscolmajor)
            {
                cblas_cgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,(int)R1,(int)R2,(int)C1,o,X1,(int)R1,X2,(int)R2,z,Y,(int)R1);
            }
            else
            {
                cblas_cgemm(CblasRowMajor,CblasNoTrans,CblasConjTrans,(int)R1,(int)R2,(int)C1,o,X1,(int)C2,X2,(int)C2,z,Y,(int)R2);
            }
        }
    }

    return 0;
}


int matmul2t_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const char iscolmajor, const char tr)
{
    const size_t N = (tr) ? C1*C2 : R1*R2;
    double x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0) {}
    else if (N<1100)
    {
        const size_t N1 = R1*C1, N2 = R2*C2;
        if (tr)
        {
            if (iscolmajor)
            {
                for (size_t c2=0; c2<C2; ++c2)
                {
                    for (size_t c1=0; c1<C1; ++c1, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t r2=0; r2<R2; ++r2, ++X1, ++X2)
                        {
                            x1r = *X1; x1i = -*++X1;
                            x2r = *X2; x2i = *++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (c1<C1-1) { X2 -= 2*R2; }
                    }
                    if (c2<C2-1) { X1 -= 2*N1; }
                }
            }
            else
            {
                for (size_t c1=0; c1<C1; ++c1)
                {
                    for (size_t c2=0; c2<C2; ++c2, X1-=2*N1, X2-=2*N2-2, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t r1=0; r1<R1; ++r1, X1+=2*C1-1, X2+=2*C2-1)
                        {
                            x1r = *X1; x1i = -*++X1;
                            x2r = *X2; x2i = *++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                    }
                    if (c1<C1-1) { X1 += 2; X2 -= 2*C2; }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (size_t r2=0; r2<R2; ++r2)
                {
                    for (size_t r1=0; r1<R1; ++r1, X1-=2*N1-2, X2-=2*N2, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t c2=0; c2<C2; ++c2, X1+=2*R1-1, X2+=2*R2-1)
                        {
                            x1r = *X1; x1i = *++X1;
                            x2r = *X2; x2i = -*++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                    }
                    if (r2<R2-1) { X1 -= 2*R1; X2 += 2; }
                }
            }
            else
            {
                for (size_t r1=0; r1<R1; ++r1)
                {
                    for (size_t r2=0; r2<R2; ++r2, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t c1=0; c1<C1; ++c1, ++X1, ++X2)
                        {
                            x1r = *X1; x1i = *++X1;
                            x2r = *X2; x2i = -*++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (r2<R2-1) { X1 -= 2*C1; }
                    }
                    if (r1<R1-1) { X2 -= 2*N2; }
                }
            }
        }
    }
    else
    {
        const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
        if (tr)
        {
            if (iscolmajor)
            {
                cblas_zgemm(CblasColMajor,CblasConjTrans,CblasNoTrans,(int)C1,(int)C2,(int)R1,o,X1,(int)R1,X2,(int)R1,z,Y,(int)C1);
            }
            else
            {
                cblas_zgemm(CblasRowMajor,CblasConjTrans,CblasNoTrans,(int)C1,(int)C2,(int)R1,o,X1,(int)C1,X2,(int)C2,z,Y,(int)C2);
            }
        }
        else
        {
            if (iscolmajor)
            {
                cblas_zgemm(CblasColMajor,CblasNoTrans,CblasConjTrans,(int)R1,(int)R2,(int)C1,o,X1,(int)R1,X2,(int)R2,z,Y,(int)R1);
            }
            else
            {
                cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasConjTrans,(int)R1,(int)R2,(int)C1,o,X1,(int)C2,X2,(int)C2,z,Y,(int)R2);
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
