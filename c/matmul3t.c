//Linear algebra function.
//Matrix multiplication of 3 input matrices X1, X2 and X3':
//Y = X1*X2*X3'.
//Multiplication order is chosen optimally:
//Y = (X1*X2)*X3', or: Y = X1*(X2*X3').

//With transpose (tr) opt, Y = X1'*X2*X3.
//Y = (X1'*X2)*X3, or: Y = X1'*(X2*X3).

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mm2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor);
int mm2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor);
int mm2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor);
int mm2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor);

int mm2t_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor, const int tr);
int mm2t_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor, const int tr);
int mm2t_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor, const int tr);
int mm2t_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor, const int tr);

int matmul3t_s (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t C1, const size_t R3, const size_t C3, const int iscolmajor, const int tr);
int matmul3t_d (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t C1, const size_t R3, const size_t C3, const int iscolmajor, const int tr);
int matmul3t_c (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t C1, const size_t R3, const size_t C3, const int iscolmajor, const int tr);
int matmul3t_z (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t C1, const size_t R3, const size_t C3, const int iscolmajor, const int tr);


int mm2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor)
{
    const size_t N = R1*C2;
    float sm2;

    if (N==0u) {}
    else if (N<1100u)
    {
        if (iscolmajor)
        {
            const size_t N1 = R1*C1;
            for (size_t c2=0u; c2<C2; ++c2)
            {
                for (size_t r1=0u; r1<R1; ++r1, X1-=N1-1u, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t r2=0u; r2<R2; ++r2, X1+=R1, ++X2)
                    {
                        sm2 = fmaf(*X1,*X2,sm2);
                    }
                    *Y = sm2;
                    if (r1<R1-1u) { X2 -= C1; }
                }
                if (c2<C2-1u) { X1 -= R1; }
            }
        }
        else
        {
            const size_t N2 = C1*C2;
            for (size_t r1=0u; r1<R1; ++r1)
            {
                for (size_t c2=0u; c2<C2; ++c2, X2-=N2-1u, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t c1=0u; c1<C1; ++c1, ++X1, X2+=C2)
                    {
                        sm2 = fmaf(*X1,*X2,sm2);
                    }
                    *Y = sm2;
                    if (c2<C2-1u) { X1 -= C1; }
                }
                if (r1<R1-1u) { X2 -= C2; }
            }
        }
    }
    else
    {
        if (iscolmajor)
        {
            cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,1.0f,X1,(int)R1,X2,(int)C1,0.0f,Y,(int)R1);
        }
        else
        {
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,1.0f,X1,(int)C1,X2,(int)C2,0.0f,Y,(int)C2);
        }
    }

    return 0;
}


int mm2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor)
{
    const size_t N = R1*C2;
    double sm2;

    if (N==0u) {}
    else if (N<1100u)
    {
        if (iscolmajor)
        {
            const size_t N1 = R1*C1;
            for (size_t c2=0u; c2<C2; ++c2)
            {
                for (size_t r1=0u; r1<R1; ++r1, X1-=N1-1u, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t r2=0u; r2<R2; ++r2, X1+=R1, ++X2)
                    {
                        sm2 = fma(*X1,*X2,sm2);
                    }
                    *Y = sm2;
                    if (r1<R1-1u) { X2 -= C1; }
                }
                if (c2<C2-1u) { X1 -= R1; }
            }
        }
        else
        {
            const size_t N2 = C1*C2;
            for (size_t r1=0u; r1<R1; ++r1)
            {
                for (size_t c2=0u; c2<C2; ++c2, X2-=N2-1u, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t c1=0u; c1<C1; ++c1, ++X1, X2+=C2)
                    {
                        sm2 = fma(*X1,*X2,sm2);
                    }
                    *Y = sm2;
                    if (c2<C2-1u) { X1 -= C1; }
                }
                if (r1<R1-1u) { X2 -= C2; }
            }
        }
    }
    else
    {
        if (iscolmajor)
        {
            cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,1.0,X1,(int)R1,X2,(int)C1,0.0,Y,(int)R1);
        }
        else
        {
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,1.0,X1,(int)C1,X2,(int)C2,0.0,Y,(int)C2);
        }
    }

    return 0;
}


int mm2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor)
{
    const size_t N = R1*C2;
    float x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0u) {}
    else if (N<1100u)
    {
        if (iscolmajor)
        {
            const size_t N1 = R1*C1;
            for (size_t c2=0u; c2<C2; ++c2)
            {
                for (size_t r1=0u; r1<R1; ++r1, X1-=2u*N1-2u, ++Y)
                {
                    sm2r = sm2i = 0.0f;
                    for (size_t r2=0u; r2<R2; ++r2, X1+=2u*R1-1u, ++X2)
                    {
                        x1r = *X1; x1i = *++X1;
                        x2r = *X2; x2i = *++X2;
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                    if (r1<R1-1u) { X2 -= 2u*C1; }
                }
                if (c2<C2-1u) { X1 -= 2u*R1; }
            }
        }
        else
        {
            const size_t N2 = C1*C2;
            for (size_t r1=0u; r1<R1; ++r1)
            {
                for (size_t c2=0u; c2<C2; ++c2, X2-=2u*N2-2u, ++Y)
                {
                    sm2r = sm2i = 0.0f;
                    for (size_t c1=0u; c1<C1; ++c1, ++X1, X2+=2u*C2-1u)
                    {
                        x1r = *X1; x1i = *++X1;
                        x2r = *X2; x2i = *++X2;
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                    if (c2<C2-1u) { X1 -= 2u*C1; }
                }
                if (r1<R1-1u) { X2 -= 2u*C2; }
            }
        }
    }
    else
    {
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
        if (iscolmajor)
        {
            cblas_cgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,o,X1,(int)R1,X2,(int)C1,z,Y,(int)R1);
        }
        else
        {
            cblas_cgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,o,X1,(int)C1,X2,(int)C2,z,Y,(int)C2);
        }
    }

    return 0;
}


int mm2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor)
{
    const size_t N = R1*C2;
    double x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0u) {}
    else if (N<1100u)
    {
        if (iscolmajor)
        {
            const size_t N1 = R1*C1;
            for (size_t c2=0u; c2<C2; ++c2)
            {
                for (size_t r1=0u; r1<R1; ++r1, X1-=2u*N1-2u, ++Y)
                {
                    sm2r = sm2i = 0.0;
                    for (size_t r2=0u; r2<R2; ++r2, X1+=2u*R1-1u, ++X2)
                    {
                        x1r = *X1; x1i = *++X1;
                        x2r = *X2; x2i = *++X2;
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                    if (r1<R1-1u) { X2 -= 2u*C1; }
                }
                if (c2<C2-1u) { X1 -= 2u*R1; }
            }
        }
        else
        {
            const size_t N2 = C1*C2;
            for (size_t r1=0u; r1<R1; ++r1)
            {
                for (size_t c2=0u; c2<C2; ++c2, X2-=2u*N2-2u, ++Y)
                {
                    sm2r = sm2i = 0.0;
                    for (size_t c1=0u; c1<C1; ++c1, ++X1, X2+=2u*C2-1u)
                    {
                        x1r = *X1; x1i = *++X1;
                        x2r = *X2; x2i = *++X2;
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                    if (c2<C2-1u) { X1 -= 2u*C1; }
                }
                if (r1<R1-1u) { X2 -= 2u*C2; }
            }
        }
    }
    else
    {
        const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
        if (iscolmajor)
        {
            cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,o,X1,(int)R1,X2,(int)C1,z,Y,(int)R1);
        }
        else
        {
            cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,o,X1,(int)C1,X2,(int)C2,z,Y,(int)C2);
        }
    }

    return 0;
}


int mm2t_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor, const int tr)
{
    const size_t N = (tr) ? C1*C2 : R1*R2;
    float sm2;

    if (N==0u) {}
    else if (N<1100u)
    {
        const size_t N1 = R1*C1, N2 = R2*C2;
        if (tr)
        {
            if (iscolmajor)
            {
                for (size_t c2=0u; c2<C2; ++c2)
                {
                    for (size_t c1=0u; c1<C1; ++c1, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t r2=0u; r2<R2; ++r2, ++X1, ++X2) { sm2 = fmaf(*X1,*X2,sm2); }
                        *Y = sm2;
                        if (c1<C1-1u) { X2 -= R2; }
                    }
                    if (c2<C2-1u) { X1 -= N1; }
                }
            }
            else
            {
                for (size_t c1=0u; c1<C1; ++c1)
                {
                    for (size_t c2=0u; c2<C2; ++c2, X1-=N1, X2-=N2-1u, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t r1=0u; r1<R1; ++r1, X1+=C1, X2+=C2) { sm2 = fmaf(*X1,*X2,sm2); }
                        *Y = sm2;
                    }
                    if (c1<C1-1u) { ++X1; X2 -= C2; }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (size_t r2=0u; r2<R2; ++r2)
                {
                    for (size_t r1=0u; r1<R1; ++r1, X1-=N1-1u, X2-=N2, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t c2=0u; c2<C2; ++c2, X1+=R1, X2+=R2) { sm2 = fmaf(*X1,*X2,sm2); }
                        *Y = sm2;
                    }
                    if (r2<R2-1u) { X1 -= R1; ++X2; }
                }
            }
            else
            {
                for (size_t r1=0u; r1<R1; ++r1)
                {
                    for (size_t r2=0u; r2<R2; ++r2, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t c1=0u; c1<C1; ++c1, ++X1, ++X2) { sm2 = fmaf(*X1,*X2,sm2); }
                        *Y = sm2;
                        if (r2<R2-1u) { X1 -= C1; }
                    }
                    if (r1<R1-1u) { X2 -= N2; }
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


int mm2t_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor, const int tr)
{
    const size_t N = (tr) ? C1*C2 : R1*R2;
    double sm2;

    if (N==0u) {}
    else if (N<1100u)
    {
        const size_t N1 = R1*C1, N2 = R2*C2;
        if (tr)
        {
            if (iscolmajor)
            {
                for (size_t c2=0u; c2<C2; ++c2)
                {
                    for (size_t c1=0u; c1<C1; ++c1, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t r2=0u; r2<R2; ++r2, ++X1, ++X2) { sm2 = fma(*X1,*X2,sm2); }
                        *Y = sm2;
                        if (c1<C1-1u) { X2 -= R2; }
                    }
                    if (c2<C2-1u) { X1 -= N1; }
                }
            }
            else
            {
                for (size_t c1=0u; c1<C1; ++c1)
                {
                    for (size_t c2=0u; c2<C2; ++c2, X1-=N1, X2-=N2-1u, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t r1=0u; r1<R1; ++r1, X1+=C1, X2+=C2) { sm2 = fma(*X1,*X2,sm2); }
                        *Y = sm2;
                    }
                    if (c1<C1-1u) { ++X1; X2 -= C2; }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (size_t r2=0u; r2<R2; ++r2)
                {
                    for (size_t r1=0u; r1<R1; ++r1, X1-=N1-1u, X2-=N2, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t c2=0u; c2<C2; ++c2, X1+=R1, X2+=R2) { sm2 = fma(*X1,*X2,sm2); }
                        *Y = sm2;
                    }
                    if (r2<R2-1u) { X1 -= R1; ++X2; }
                }
            }
            else
            {
                for (size_t r1=0u; r1<R1; ++r1)
                {
                    for (size_t r2=0u; r2<R2; ++r2, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t c1=0u; c1<C1; ++c1, ++X1, ++X2) { sm2 = fma(*X1,*X2,sm2); }
                        *Y = sm2;
                        if (r2<R2-1u) { X1 -= C1; }
                    }
                    if (r1<R1-1u) { X2 -= N2; }
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


int mm2t_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor, const int tr)
{
    const size_t N = (tr) ? C1*C2 : R1*R2;
    float x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0u) {}
    else if (N<1100u)
    {
        const size_t N1 = R1*C1, N2 = R2*C2;
        if (tr)
        {
            if (iscolmajor)
            {
                for (size_t c2=0u; c2<C2; ++c2)
                {
                    for (size_t c1=0u; c1<C1; ++c1, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t r2=0u; r2<R2; ++r2, ++X1, ++X2)
                        {
                            x1r = *X1; x1i = -*++X1;
                            x2r = *X2; x2i = *++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (c1<C1-1u) { X2 -= 2u*R2; }
                    }
                    if (c2<C2-1u) { X1 -= 2u*N1; }
                }
            }
            else
            {
                for (size_t c1=0u; c1<C1; ++c1)
                {
                    for (size_t c2=0u; c2<C2; ++c2, X1-=2u*N1, X2-=2u*N2-2u, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t r1=0u; r1<R1; ++r1, X1+=2u*C1-1u, X2+=2u*C2-1u)
                        {
                            x1r = *X1; x1i = -*++X1;
                            x2r = *X2; x2i = *++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                    }
                    if (c1<C1-1u) { X1 += 2; X2 -= 2u*C2; }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (size_t r2=0u; r2<R2; ++r2)
                {
                    for (size_t r1=0u; r1<R1; ++r1, X1-=2u*N1-2u, X2-=2u*N2, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t c2=0u; c2<C2; ++c2, X1+=2u*R1-1u, X2+=2u*R2-1u)
                        {
                            x1r = *X1; x1i = *++X1;
                            x2r = *X2; x2i = -*++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                    }
                    if (r2<R2-1u) { X1 -= 2u*R1; X2 += 2; }
                }
            }
            else
            {
                for (size_t r1=0u; r1<R1; ++r1)
                {
                    for (size_t r2=0u; r2<R2; ++r2, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t c1=0u; c1<C1; ++c1, ++X1, ++X2)
                        {
                            x1r = *X1; x1i = *++X1;
                            x2r = *X2; x2i = -*++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (r2<R2-1u) { X1 -= 2u*C1; }
                    }
                    if (r1<R1-1u) { X2 -= 2u*N2; }
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


int mm2t_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor, const int tr)
{
    const size_t N = (tr) ? C1*C2 : R1*R2;
    double x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0u) {}
    else if (N<1100u)
    {
        const size_t N1 = R1*C1, N2 = R2*C2;
        if (tr)
        {
            if (iscolmajor)
            {
                for (size_t c2=0u; c2<C2; ++c2)
                {
                    for (size_t c1=0u; c1<C1; ++c1, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t r2=0u; r2<R2; ++r2, ++X1, ++X2)
                        {
                            x1r = *X1; x1i = -*++X1;
                            x2r = *X2; x2i = *++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (c1<C1-1u) { X2 -= 2u*R2; }
                    }
                    if (c2<C2-1u) { X1 -= 2u*N1; }
                }
            }
            else
            {
                for (size_t c1=0u; c1<C1; ++c1)
                {
                    for (size_t c2=0u; c2<C2; ++c2, X1-=2u*N1, X2-=2u*N2-2u, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t r1=0u; r1<R1; ++r1, X1+=2u*C1-1u, X2+=2u*C2-1u)
                        {
                            x1r = *X1; x1i = -*++X1;
                            x2r = *X2; x2i = *++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                    }
                    if (c1<C1-1u) { X1 += 2; X2 -= 2u*C2; }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (size_t r2=0u; r2<R2; ++r2)
                {
                    for (size_t r1=0u; r1<R1; ++r1, X1-=2u*N1-2u, X2-=2u*N2, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t c2=0u; c2<C2; ++c2, X1+=2u*R1-1u, X2+=2u*R2-1u)
                        {
                            x1r = *X1; x1i = *++X1;
                            x2r = *X2; x2i = -*++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                    }
                    if (r2<R2-1u) { X1 -= 2u*R1; X2 += 2; }
                }
            }
            else
            {
                for (size_t r1=0u; r1<R1; ++r1)
                {
                    for (size_t r2=0u; r2<R2; ++r2, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t c1=0u; c1<C1; ++c1, ++X1, ++X2)
                        {
                            x1r = *X1; x1i = *++X1;
                            x2r = *X2; x2i = -*++X2;
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (r2<R2-1u) { X1 -= 2u*C1; }
                    }
                    if (r1<R1-1u) { X2 -= 2u*N2; }
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


int matmul3t_s (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t C1, const size_t R3, const size_t C3, const int iscolmajor, const int tr)
{
    if (tr)
    {
        const size_t R2 = R1, C2 = R3;
        if (C1*C3==0u) {}
        else if (R1*C1*(C2+C3)<R2*C3*(C1+R3))
        {
            float *X12;
            if (!(X12=(float *)malloc(C1*C2*sizeof(float)))) { fprintf(stderr,"error in matmul3t_s: problem with malloc. "); perror("malloc"); return 1; }
            mm2t_s(X12,X1,X2,R1,C1,R2,C2,iscolmajor,tr);
            mm2_s(Y,X12,X3,C1,C2,R3,C3,iscolmajor);
            free(X12);
        }
        else
        {
            float *X23;
            if (!(X23=(float *)malloc(R2*C3*sizeof(float)))) { fprintf(stderr,"error in matmul3t_s: problem with malloc. "); perror("malloc"); return 1; }
            mm2_s(X23,X2,X3,R2,C2,R3,C3,iscolmajor);
            mm2t_s(Y,X1,X23,R1,C1,R2,C3,iscolmajor,tr);
            free(X23);
        }
    }
    else
    {
        const size_t R2 = C1, C2 = C3;
        if (R1*R3==0u) {}
        else if (R1*C2*(C1+R3)<C1*R3*(R1+C3))
        {
            float *X12;
            if (!(X12=(float *)malloc(R1*C2*sizeof(float)))) { fprintf(stderr,"error in matmul3t_s: problem with malloc. "); perror("malloc"); return 1; }
            mm2_s(X12,X1,X2,R1,C1,R2,C2,iscolmajor);
            mm2t_s(Y,X12,X3,R1,C2,R3,C3,iscolmajor,tr);
            free(X12);
        }
        else
        {
            float *X23;
            if (!(X23=(float *)malloc(R2*R3*sizeof(float)))) { fprintf(stderr,"error in matmul3t_s: problem with malloc. "); perror("malloc"); return 1; }
            mm2t_s(X23,X2,X3,R2,C2,R3,C3,iscolmajor,tr);
            mm2_s(Y,X1,X23,R1,C1,R2,R3,iscolmajor);
            free(X23);
        }
    }

    return 0;
}


int matmul3t_d (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t C1, const size_t R3, const size_t C3, const int iscolmajor, const int tr)
{
    if (tr)
    {
        const size_t R2 = R1, C2 = R3;
        if (C1*C3==0u) {}
        else if (R1*C1*(C2+C3)<R2*C3*(C1+R3))
        {
            double *X12;
            if (!(X12=(double *)malloc(C1*C2*sizeof(double)))) { fprintf(stderr,"error in matmul3t_d: problem with malloc. "); perror("malloc"); return 1; }
            mm2t_d(X12,X1,X2,R1,C1,R2,C2,iscolmajor,tr);
            mm2_d(Y,X12,X3,C1,C2,R3,C3,iscolmajor);
            free(X12);
        }
        else
        {
            double *X23;
            if (!(X23=(double *)malloc(R2*C3*sizeof(double)))) { fprintf(stderr,"error in matmul3t_d: problem with malloc. "); perror("malloc"); return 1; }
            mm2_d(X23,X2,X3,R2,C2,R3,C3,iscolmajor);
            mm2t_d(Y,X1,X23,R1,C1,R2,C3,iscolmajor,tr);
            free(X23);
        }
    }
    else
    {
        const size_t R2 = C1, C2 = C3;
        if (R1*R3==0u) {}
        else if (R1*C2*(C1+R3)<C1*R3*(R1+C3))
        {
            double *X12;
            if (!(X12=(double *)malloc(R1*C2*sizeof(double)))) { fprintf(stderr,"error in matmul3t_d: problem with malloc. "); perror("malloc"); return 1; }
            mm2_d(X12,X1,X2,R1,C1,R2,C2,iscolmajor);
            mm2t_d(Y,X12,X3,R1,C2,R3,C3,iscolmajor,tr);
            free(X12);
        }
        else
        {
            double *X23;
            if (!(X23=(double *)malloc(R2*R3*sizeof(double)))) { fprintf(stderr,"error in matmul3t_d: problem with malloc. "); perror("malloc"); return 1; }
            mm2t_d(X23,X2,X3,R2,C2,R3,C3,iscolmajor,tr);
            mm2_d(Y,X1,X23,R1,C1,R2,R3,iscolmajor);
            free(X23);
        }
    }

    return 0;
}


int matmul3t_c (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t C1, const size_t R3, const size_t C3, const int iscolmajor, const int tr)
{
    if (tr)
    {
        const size_t R2 = R1, C2 = R3;
        if (C1*C3==0u) {}
        else if (R1*C1*(C2+C3)<R2*C3*(C1+R3))
        {
            float *X12;
            if (!(X12=(float *)malloc(2u*C1*C2*sizeof(float)))) { fprintf(stderr,"error in matmul3t_c: problem with malloc. "); perror("malloc"); return 1; }
            mm2t_c(X12,X1,X2,R1,C1,R2,C2,iscolmajor,tr);
            mm2_c(Y,X12,X3,C1,C2,R3,C3,iscolmajor);
            free(X12);
        }
        else
        {
            float *X23;
            if (!(X23=(float *)malloc(2u*R2*C3*sizeof(float)))) { fprintf(stderr,"error in matmul3t_c: problem with malloc. "); perror("malloc"); return 1; }
            mm2_c(X23,X2,X3,R2,C2,R3,C3,iscolmajor);
            mm2t_c(Y,X1,X23,R1,C1,R2,C3,iscolmajor,tr);
            free(X23);
        }
    }
    else
    {
        const size_t R2 = C1, C2 = C3;
        if (R1*R3==0u) {}
        else if (R1*C2*(C1+R3)<C1*R3*(R1+C3))
        {
            float *X12;
            if (!(X12=(float *)malloc(2u*R1*C2*sizeof(float)))) { fprintf(stderr,"error in matmul3t_c: problem with malloc. "); perror("malloc"); return 1; }
            mm2_c(X12,X1,X2,R1,C1,R2,C2,iscolmajor);
            mm2t_c(Y,X12,X3,R1,C2,R3,C3,iscolmajor,tr);
            free(X12);
        }
        else
        {
            float *X23;
            if (!(X23=(float *)malloc(2u*R2*R3*sizeof(float)))) { fprintf(stderr,"error in matmul3t_c: problem with malloc. "); perror("malloc"); return 1; }
            mm2t_c(X23,X2,X3,R2,C2,R3,C3,iscolmajor,tr);
            mm2_c(Y,X1,X23,R1,C1,R2,R3,iscolmajor);
            free(X23);
        }
    }

    return 0;
}


int matmul3t_z (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t C1, const size_t R3, const size_t C3, const int iscolmajor, const int tr)
{
    if (tr)
    {
        const size_t R2 = R1, C2 = R3;
        if (C1*C3==0u) {}
        else if (R1*C1*(C2+C3)<R2*C3*(C1+R3))
        {
            double *X12;
            if (!(X12=(double *)malloc(2u*C1*C2*sizeof(double)))) { fprintf(stderr,"error in matmul3t_z: problem with malloc. "); perror("malloc"); return 1; }
            mm2t_z(X12,X1,X2,R1,C1,R2,C2,iscolmajor,tr);
            mm2_z(Y,X12,X3,C1,C2,R3,C3,iscolmajor);
            free(X12);
        }
        else
        {
            double *X23;
            if (!(X23=(double *)malloc(2u*R2*C3*sizeof(double)))) { fprintf(stderr,"error in matmul3t_z: problem with malloc. "); perror("malloc"); return 1; }
            mm2_z(X23,X2,X3,R2,C2,R3,C3,iscolmajor);
            mm2t_z(Y,X1,X23,R1,C1,R2,C3,iscolmajor,tr);
            free(X23);
        }
    }
    else
    {
        const size_t R2 = C1, C2 = C3;
        if (R1*R3==0u) {}
        else if (R1*C2*(C1+R3)<C1*R3*(R1+C3))
        {
            double *X12;
            if (!(X12=(double *)malloc(2u*R1*C2*sizeof(double)))) { fprintf(stderr,"error in matmul3t_z: problem with malloc. "); perror("malloc"); return 1; }
            mm2_z(X12,X1,X2,R1,C1,R2,C2,iscolmajor);
            mm2t_z(Y,X12,X3,R1,C2,R3,C3,iscolmajor,tr);
            free(X12);
        }
        else
        {
            double *X23;
            if (!(X23=(double *)malloc(2u*R2*R3*sizeof(double)))) { fprintf(stderr,"error in matmul3t_z: problem with malloc. "); perror("malloc"); return 1; }
            mm2t_z(X23,X2,X3,R2,C2,R3,C3,iscolmajor,tr);
            mm2_z(Y,X1,X23,R1,C1,R2,R3,iscolmajor);
            free(X23);
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
