//Linear algebra function.
//Matrix multiplication of 3 input matrices X1, X2 and X3:
//Y = X1*X2*X3.
//Multiplication order is chosen optimally:
//Y = (X1*X2)*X3
//or:
//Y = X1*(X2*X3)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mm2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor);
int mm2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor);
int mm2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor);
int mm2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor);

int matmul3_s (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t C1, const size_t C2, const size_t C3, const char iscolmajor);
int matmul3_d (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t C1, const size_t C2, const size_t C3, const char iscolmajor);
int matmul3_c (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t C1, const size_t C2, const size_t C3, const char iscolmajor);
int matmul3_z (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t C1, const size_t C2, const size_t C3, const char iscolmajor);


int mm2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor)
{
    const size_t N = R1*C2;
    float sm2;

    if (N==0) {}
    else if (N<1100)
    {
        if (iscolmajor)
        {
            const size_t N1 = R1*C1;
            for (size_t c2=0; c2<C2; ++c2)
            {
                for (size_t r1=0; r1<R1; ++r1, X1-=N1-1, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t r2=0; r2<C1; ++r2, X1+=R1, ++X2)
                    {
                        sm2 = fmaf(*X1,*X2,sm2);
                    }
                    *Y = sm2;
                    if (r1<R1-1) { X2 -= C1; }
                }
                if (c2<C2-1) { X1 -= R1; }
            }
        }
        else
        {
            const size_t N2 = C1*C2;
            for (size_t r1=0; r1<R1; ++r1)
            {
                for (size_t c2=0; c2<C2; ++c2, X2-=N2-1, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t c1=0; c1<C1; ++c1, ++X1, X2+=C2)
                    {
                        sm2 = fmaf(*X1,*X2,sm2);
                    }
                    *Y = sm2;
                    if (c2<C2-1) { X1 -= C1; }
                }
                if (r1<R1-1) { X2 -= C2; }
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


int mm2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor)
{
    const size_t N = R1*C2;
    double sm2;

    if (N==0) {}
    else if (N<1100)
    {
        if (iscolmajor)
        {
            const size_t N1 = R1*C1;
            for (size_t c2=0; c2<C2; ++c2)
            {
                for (size_t r1=0; r1<R1; ++r1, X1-=N1-1, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t r2=0; r2<C1; ++r2, X1+=R1, ++X2)
                    {
                        sm2 = fma(*X1,*X2,sm2);
                    }
                    *Y = sm2;
                    if (r1<R1-1) { X2 -= C1; }
                }
                if (c2<C2-1) { X1 -= R1; }
            }
        }
        else
        {
            const size_t N2 = C1*C2;
            for (size_t r1=0; r1<R1; ++r1)
            {
                for (size_t c2=0; c2<C2; ++c2, X2-=N2-1, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t c1=0; c1<C1; ++c1, ++X1, X2+=C2)
                    {
                        sm2 = fma(*X1,*X2,sm2);
                    }
                    *Y = sm2;
                    if (c2<C2-1) { X1 -= C1; }
                }
                if (r1<R1-1) { X2 -= C2; }
            }
        }
    }
    else
    {
        if (iscolmajor)
        {
            cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,1.0,X1,(int)R1,X2,(int)C1,1.0,Y,(int)R1);
        }
        else
        {
            cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,1.0,X1,(int)C1,X2,(int)C2,1.0,Y,(int)C2);
        }
    }

    return 0;
}


int mm2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor)
{
    const size_t N = R1*C2;
    float x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0) {}
    else if (N<1100)
    {
        if (iscolmajor)
        {
            const size_t N1 = R1*C1;
            for (size_t c2=0; c2<C2; ++c2)
            {
                for (size_t r1=0; r1<R1; ++r1, X1-=2*N1-2, ++Y)
                {
                    sm2r = sm2i = 0.0f;
                    for (size_t r2=0; r2<C1; ++r2, X1+=2*R1-1, ++X2)
                    {
                        x1r = *X1; x1i = *++X1;
                        x2r = *X2; x2i = *++X2;
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                    if (r1<R1-1) { X2 -= 2*C1; }
                }
                if (c2<C2-1) { X1 -= 2*R1; }
            }
        }
        else
        {
            const size_t N2 = C1*C2;
            for (size_t r1=0; r1<R1; ++r1)
            {
                for (size_t c2=0; c2<C2; ++c2, X2-=2*N2-2, ++Y)
                {
                    sm2r = sm2i = 0.0f;
                    for (size_t c1=0; c1<C1; ++c1, ++X1, X2+=2*C2-1)
                    {
                        x1r = *X1; x1i = *++X1;
                        x2r = *X2; x2i = *++X2;
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                    if (c2<C2-1) { X1 -= 2*C1; }
                }
                if (r1<R1-1) { X2 -= 2*C2; }
            }
        }
    }
    else
    {
        const float o[2] = {1.0f,0.0f};
        if (iscolmajor)
        {
            cblas_cgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,o,X1,(int)R1,X2,(int)C1,o,Y,(int)R1);
        }
        else
        {
            cblas_cgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,o,X1,(int)C1,X2,(int)C2,o,Y,(int)C2);
        }
    }

    return 0;
}


int mm2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t C2, const char iscolmajor)
{
    const size_t N = R1*C2;
    double x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0) {}
    else if (N<1100)
    {
        if (iscolmajor)
        {
            const size_t N1 = R1*C1;
            for (size_t c2=0; c2<C2; ++c2)
            {
                for (size_t r1=0; r1<R1; ++r1, X1-=2*N1-2, ++Y)
                {
                    sm2r = sm2i = 0.0;
                    for (size_t r2=0; r2<C1; ++r2, X1+=2*R1-1, ++X2)
                    {
                        x1r = *X1; x1i = *++X1;
                        x2r = *X2; x2i = *++X2;
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                    if (r1<R1-1) { X2 -= 2*C1; }
                }
                if (c2<C2-1) { X1 -= 2*R1; }
            }
        }
        else
        {
            const size_t N2 = C1*C2;
            for (size_t r1=0; r1<R1; ++r1)
            {
                for (size_t c2=0; c2<C2; ++c2, X2-=2*N2-2, ++Y)
                {
                    sm2r = sm2i = 0.0;
                    for (size_t c1=0; c1<C1; ++c1, ++X1, X2+=2*C2-1)
                    {
                        x1r = *X1; x1i = *++X1;
                        x2r = *X2; x2i = *++X2;
                        sm2r += x1r*x2r - x1i*x2i;
                        sm2i += x1r*x2i + x1i*x2r;
                    }
                    *Y = sm2r; *++Y = sm2i;
                    if (c2<C2-1) { X1 -= 2*C1; }
                }
                if (r1<R1-1) { X2 -= 2*C2; }
            }
        }
    }
    else
    {
        const double o[2] = {1.0,0.0};
        if (iscolmajor)
        {
            cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,o,X1,(int)R1,X2,(int)C1,o,Y,(int)R1);
        }
        else
        {
            cblas_zgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,(int)R1,(int)C2,(int)C1,o,X1,(int)C1,X2,(int)C2,o,Y,(int)C2);
        }
    }

    return 0;
}


int matmul3_s (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t C1, const size_t C2, const size_t C3, const char iscolmajor)
{
    if (R1*C3==0) {}
    else if (R1*C2*(C1+C3)<C1*C3*(R1+C2))
    {
        float *X12;
        if (!(X12=(float *)malloc(R1*C2*sizeof(float)))) { fprintf(stderr,"error in matmul3_s: problem with malloc. "); perror("malloc"); return 1; }
        mm2_s(X12,X1,X2,R1,C1,C2,iscolmajor);
        mm2_s(Y,X12,X3,R1,C2,C3,iscolmajor);
        free(X12);
    }
    else
    {
        float *X23;
        if (!(X23=(float *)malloc(C1*C3*sizeof(float)))) { fprintf(stderr,"error in matmul3_s: problem with malloc. "); perror("malloc"); return 1; }
        mm2_s(X23,X2,X3,C1,C2,C3,iscolmajor);
        mm2_s(Y,X1,X23,R1,C1,C3,iscolmajor);
        free(X23);
    }

    return 0;
}


int matmul3_d (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t C1, const size_t C2, const size_t C3, const char iscolmajor)
{
    if (R1*C3==0) {}
    else if (R1*C2*(C1+C3)<C1*C3*(R1+C2))
    {
        double *X12;
        if (!(X12=(double *)malloc(R1*C2*sizeof(double)))) { fprintf(stderr,"error in matmul3_d: problem with malloc. "); perror("malloc"); return 1; }
        mm2_d(X12,X1,X2,R1,C1,C2,iscolmajor);
        mm2_d(Y,X12,X3,R1,C2,C3,iscolmajor);
        free(X12);
    }
    else
    {
        double *X23;
        if (!(X23=(double *)malloc(C1*C3*sizeof(double)))) { fprintf(stderr,"error in matmul3_d: problem with malloc. "); perror("malloc"); return 1; }
        mm2_d(X23,X2,X3,C1,C2,C3,iscolmajor);
        mm2_d(Y,X1,X23,R1,C1,C3,iscolmajor);
        free(X23);
    }

    return 0;
}


int matmul3_c (float *Y, const float *X1, const float *X2, const float *X3, const size_t R1, const size_t C1, const size_t C2, const size_t C3, const char iscolmajor)
{
    if (R1*C3==0) {}
    else if (R1*C2*(C1+C3)<C1*C3*(R1+C2))
    {
        float *X12;
        if (!(X12=(float *)malloc(2*R1*C2*sizeof(float)))) { fprintf(stderr,"error in matmul3_c: problem with malloc. "); perror("malloc"); return 1; }
        mm2_c(X12,X1,X2,R1,C1,C2,iscolmajor);
        mm2_c(Y,X12,X3,R1,C2,C3,iscolmajor);
        free(X12);
    }
    else
    {
        float *X23;
        if (!(X23=(float *)malloc(2*C1*C3*sizeof(float)))) { fprintf(stderr,"error in matmul3_c: problem with malloc. "); perror("malloc"); return 1; }
        mm2_c(X23,X2,X3,C1,C2,C3,iscolmajor);
        mm2_c(Y,X1,X23,R1,C1,C3,iscolmajor);
        free(X23);
    }

    return 0;
}


int matmul3_z (double *Y, const double *X1, const double *X2, const double *X3, const size_t R1, const size_t C1, const size_t C2, const size_t C3, const char iscolmajor)
{
    if (R1*C3==0) {}
    else if (R1*C2*(C1+C3)<C1*C3*(R1+C2))
    {
        double *X12;
        if (!(X12=(double *)malloc(2*R1*C2*sizeof(double)))) { fprintf(stderr,"error in matmul3_z: problem with malloc. "); perror("malloc"); return 1; }
        mm2_z(X12,X1,X2,R1,C1,C2,iscolmajor);
        mm2_z(Y,X12,X3,R1,C2,C3,iscolmajor);
        free(X12);
    }
    else
    {
        double *X23;
        if (!(X23=(double *)malloc(2*C1*C3*sizeof(double)))) { fprintf(stderr,"error in matmul3_z: problem with malloc. "); perror("malloc"); return 1; }
        mm2_z(X23,X2,X3,C1,C2,C3,iscolmajor);
        mm2_z(Y,X1,X23,R1,C1,C3,iscolmajor);
        free(X23);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
