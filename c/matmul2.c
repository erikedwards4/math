//Linear algebra function.
//Matrix multiplication of 2 input matrices X1 and X2.

//Note that the N<1100u applies to using fmaf and fma,
//but it is usually faster to do sm += x1 * x2.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int matmul2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor)
{
    if (C1!=R2) { fprintf(stderr,"error in matmul2_s: C1 (ncols X1) must equal R2 (nrows X2) \n"); return 1; }
    
    const size_t N = R1*C2;
    float sm2;

    if (N==0u) {}
    else if (N<1100u)
    {
        if (iscolmajor)
        {
            const size_t N1 = R1*C1;
            for (size_t c2=C2; c2>0u; --c2)
            {
                for (size_t r1=R1; r1>0u; --r1, X1-=N1-1u, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t r2=R2; r2>0u; --r2, X1+=R1, ++X2)
                    {
                        //sm2 += *X1 * *X2;
                        sm2 = fmaf(*X1,*X2,sm2);
                    }
                    *Y = sm2;
                    if (r1>1u) { X2 -= R2; }
                }
                if (c2>1u) { X1 -= R1; }
            }
        }
        else
        {
            const size_t N2 = R2*C2;
            for (size_t r1=R1; r1>0u; --r1)
            {
                for (size_t c2=C2; c2>0u; --c2, X2-=N2-1u, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t c1=C1; c1>0u; --c1, ++X1, X2+=C2)
                    {
                        //sm2 += *X1 * *X2;
                        sm2 = fmaf(*X1,*X2,sm2);
                    }
                    *Y = sm2;
                    if (c2>1u) { X1 -= C1; }
                }
                if (r1>1u) { X2 -= C2; }
            }
        }
    }
    else
    {
        if (iscolmajor)
        {
            cblas_sgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)R2,1.0f,X1,(int)R1,X2,(int)R2,0.0f,Y,(int)R1);
        }
        else
        {
            cblas_sgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,1.0f,X1,(int)C1,X2,(int)C2,0.0f,Y,(int)C2);
        }
    }

    return 0;
}


int matmul2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor)
{
    if (C1!=R2) { fprintf(stderr,"error in matmul2_d: C1 (ncols X1) must equal R2 (nrows X2) \n"); return 1; }
    
    const size_t N = R1*C2;
    double sm2;

    if (N==0u) {}
    else if (N<1100u)
    {
        if (iscolmajor)
        {
            const size_t N1 = R1*C1;
            for (size_t c2=C2; c2>0u; --c2)
            {
                for (size_t r1=R1; r1>0u; --r1, X1-=N1-1u, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t r2=R2; r2>0u; --r2, X1+=R1, ++X2)
                    {
                        sm2 = fma(*X1,*X2,sm2);
                    }
                    *Y = sm2;
                    if (r1>1u) { X2 -= R2; }
                }
                if (c2>1u) { X1 -= R1; }
            }
        }
        else
        {
            const size_t N2 = R2*C2;
            for (size_t r1=R1; r1>0u; --r1)
            {
                for (size_t c2=C2; c2>0u; --c2, X2-=N2-1u, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t c1=C1; c1>0u; --c1, ++X1, X2+=C2)
                    {
                        sm2 = fma(*X1,*X2,sm2);
                    }
                    *Y = sm2;
                    if (c2>1u) { X1 -= C1; }
                }
                if (r1>1u) { X2 -= C2; }
            }
        }
    }
    else
    {
        if (iscolmajor)
        {
            cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)R2,1.0,X1,(int)R1,X2,(int)R2,0.0,Y,(int)R1);
        }
        else
        {
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,1.0,X1,(int)C1,X2,(int)C2,0.0,Y,(int)C2);
        }
    }

    return 0;
}


int matmul2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor)
{
    if (C1!=R2) { fprintf(stderr,"error in matmul2_c: C1 (ncols X1) must equal R2 (nrows X2) \n"); return 1; }
    
    const size_t N = R1*C2;
    float x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0u) {}
    else if (N<1100u)
    {
        if (iscolmajor)
        {
            const size_t N1 = R1*C1;
            for (size_t c2=C2; c2>0u; --c2)
            {
                for (size_t r1=R1; r1>0u; --r1, X1-=2u*N1-2u, ++Y)
                {
                    sm2r = sm2i = 0.0f;
                    for (size_t r2=R2; r2>0u; --r2, X1+=2u*R1-1u, ++X2)
                    {
                        x1r = *X1; x1i = *++X1;
                        x2r = *X2; x2i = *++X2;
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                    if (r1>1u) { X2 -= 2u*R2; }
                }
                if (c2>1u) { X1 -= 2u*R1; }
            }
        }
        else
        {
            const size_t N2 = R2*C2;
            for (size_t r1=R1; r1>0u; --r1)
            {
                for (size_t c2=C2; c2>0u; --c2, X2-=2u*N2-2u, ++Y)
                {
                    sm2r = sm2i = 0.0f;
                    for (size_t c1=C1; c1>0u; --c1, ++X1, X2+=2u*C2-1u)
                    {
                        x1r = *X1; x1i = *++X1;
                        x2r = *X2; x2i = *++X2;
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                    if (c2>1u) { X1 -= 2u*C1; }
                }
                if (r1>1u) { X2 -= 2u*C2; }
            }
        }
    }
    else
    {
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
        if (iscolmajor)
        {
            cblas_cgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)R2,o,X1,(int)R1,X2,(int)R2,z,Y,(int)R1);
        }
        else
        {
            cblas_cgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,o,X1,(int)C1,X2,(int)C2,z,Y,(int)C2);
        }
    }

    return 0;
}


int matmul2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t R2, const size_t C2, const int iscolmajor)
{
    if (C1!=R2) { fprintf(stderr,"error in matmul2_z: C1 (ncols X1) must equal R2 (nrows X2) \n"); return 1; }
    
    const size_t N = R1*C2;
    double x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0u) {}
    else if (N<1100u)
    {
        if (iscolmajor)
        {
            const size_t N1 = R1*C1;
            for (size_t c2=C2; c2>0u; --c2)
            {
                for (size_t r1=R1; r1>0u; --r1, X1-=2u*N1-2u, ++Y)
                {
                    sm2r = sm2i = 0.0;
                    for (size_t r2=R2; r2>0u; --r2, X1+=2u*R1-1u, ++X2)
                    {
                        x1r = *X1; x1i = *++X1;
                        x2r = *X2; x2i = *++X2;
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                    if (r1>1u) { X2 -= 2u*R2; }
                }
                if (c2>1u) { X1 -= 2u*R1; }
            }
        }
        else
        {
            const size_t N2 = R2*C2;
            for (size_t r1=R1; r1>0u; --r1)
            {
                for (size_t c2=C2; c2>0u; --c2, X2-=2u*N2-2u, ++Y)
                {
                    sm2r = sm2i = 0.0;
                    for (size_t c1=C1; c1>0u; --c1, ++X1, X2+=2u*C2-1u)
                    {
                        x1r = *X1; x1i = *++X1;
                        x2r = *X2; x2i = *++X2;
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                    if (c2>1u) { X1 -= 2u*C1; }
                }
                if (r1>1u) { X2 -= 2u*C2; }
            }
        }
    }
    else
    {
        const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
        if (iscolmajor)
        {
            cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)R2,o,X1,(int)R1,X2,(int)R2,z,Y,(int)R1);
        }
        else
        {
            cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,o,X1,(int)C1,X2,(int)C2,z,Y,(int)C2);
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
