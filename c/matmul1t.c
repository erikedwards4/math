//Linear algebra function.
//Multiplies a square matrix by itself transposed: Y = X*X'.

//For complex case, this uses the Hermitian (conjugate) transpose.
//But could not get cblas_cgemm and cblas_zgemm working!
//So, only manual solutions here for complex case (optimal up to N=2500 anyway).

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int matmul1t_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const char tr);
int matmul1t_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const char tr);
int matmul1t_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const char tr);
int matmul1t_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const char tr);


int matmul1t_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const char tr)
{
    const size_t N = R*C;
    float sm2;

    if (N==0u) {}
    else if (N<2500u)
    {
        if (tr)
        {
            if (iscolmajor)
            {
                for (size_t c2=0u; c2<C; ++c2)
                {
                    Y += c2;
                    for (size_t c1=c2; c1<C; ++c1, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t r=0u; r<R; ++r, ++X) { sm2 = fmaf(*X,*(X+(c1-c2)*R),sm2); }
                        *Y = *(Y+(int)((c1-c2)*C+c2)-(int)c1) = sm2;
                        if (c1<C-1u) { X -= R; }
                    }
                }
            }
            else
            {
                for (size_t c1=0u; c1<C; ++c1, ++X)
                {
                    Y += c1;
                    for (size_t c2=c1; c2<C; ++c2, X-=N, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t r=0u; r<R; ++r, X+=C) { sm2 = fmaf(*X,*(X+c2-c1),sm2); }
                        *Y = *(Y+(int)((c2-c1)*C+c1)-(int)c2) = sm2;
                    }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (size_t r2=0u; r2<R; ++r2, ++X)
                {
                    Y += r2;
                    for (size_t r1=r2; r1<R; ++r1, X-=N, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t c=0u; c<C; ++c, X+=R) { sm2 = fmaf(*X,*(X+r1-r2),sm2); }
                        *Y = *(Y+(int)((r1-r2)*R+r2)-(int)r1) = sm2;
                    }
                }
            }
            else
            {
                for (size_t r1=0u; r1<R; ++r1)
                {
                    Y += r1;
                    for (size_t r2=r1; r2<R; ++r2, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t c=0u; c<C; ++c, ++X) { sm2 = fmaf(*X,*(X+(r2-r1)*C),sm2); }
                        *Y = *(Y+(int)((r2-r1)*R+r1)-(int)r2) = sm2;
                        if (r2<R-1u) { X -= C; }
                    }
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
                cblas_ssyrk(CblasColMajor,CblasUpper,CblasTrans,(int)C,(int)R,1.0f,X,(int)R,0.0f,Y,(int)C);
                for (int c=0; c<(int)C-1; ++c)
                {
                    Y += c + 1;
                    for (int r=c+1; r<(int)C; ++r, ++Y) { *Y = *(Y+(int)C*(r-c)+c-r); }
                }
            }
            else
            {
                cblas_ssyrk(CblasRowMajor,CblasLower,CblasTrans,(int)C,(int)R,1.0f,X,(int)C,0.0f,Y,(int)C);
                for (int r=0; r<(int)C-1; ++r)
                {
                    Y += r + 1;
                    for (int c=r+1; c<(int)C; ++c, ++Y) { *Y = *(Y+(int)C*(c-r)+r-c); }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                cblas_ssyrk(CblasColMajor,CblasUpper,CblasNoTrans,(int)R,(int)C,1.0f,X,(int)R,0.0f,Y,(int)R);
                for (int c=0; c<(int)R-1; ++c)
                {
                    Y += c + 1;
                    for (int r=c+1; r<(int)R; ++r, ++Y) { *Y = *(Y+(int)R*(r-c)+c-r); }
                }
            }
            else
            {
                cblas_ssyrk(CblasRowMajor,CblasLower,CblasNoTrans,(int)R,(int)C,1.0f,X,(int)C,0.0f,Y,(int)R);
                for (int r=0; r<(int)R-1; ++r)
                {
                    Y += r + 1;
                    for (int c=r+1; c<(int)R; ++c, ++Y) { *Y = *(Y+(int)R*(c-r)+r-c); }
                }
            }
        }
    }

    return 0;
}


int matmul1t_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const char tr)
{
    const size_t N = R*C;
    double sm2;

    if (N==0u) {}
    else if (N<2500u)
    {
        if (tr)
        {
            if (iscolmajor)
            {
                for (size_t c2=0u; c2<C; ++c2)
                {
                    Y += c2;
                    for (size_t c1=c2; c1<C; ++c1, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t r=0u; r<R; ++r, ++X) { sm2 = fma(*X,*(X+(c1-c2)*R),sm2); }
                        *Y = *(Y+(int)((c1-c2)*C+c2)-(int)c1) = sm2;
                        if (c1<C-1u) { X -= R; }
                    }
                }
            }
            else
            {
                for (size_t c1=0u; c1<C; ++c1, ++X)
                {
                    Y += c1;
                    for (size_t c2=c1; c2<C; ++c2, X-=N, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t r=0u; r<R; ++r, X+=C) { sm2 = fma(*X,*(X+c2-c1),sm2); }
                        *Y = *(Y+(int)((c2-c1)*C+c1)-(int)c2) = sm2;
                    }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (size_t r2=0u; r2<R; ++r2, ++X)
                {
                    Y += r2;
                    for (size_t r1=r2; r1<R; ++r1, X-=N, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t c=0u; c<C; ++c, X+=R) { sm2 = fma(*X,*(X+r1-r2),sm2); }
                        *Y = *(Y+(int)((r1-r2)*R+r2)-(int)r1) = sm2;
                    }
                }
            }
            else
            {
                for (size_t r1=0u; r1<R; ++r1)
                {
                    Y += r1;
                    for (size_t r2=r1; r2<R; ++r2, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t c=0u; c<C; ++c, ++X) { sm2 = fma(*X,*(X+(r2-r1)*C),sm2); }
                        *Y = *(Y+(int)((r2-r1)*R+r1)-(int)r2) = sm2;
                        if (r2<R-1u) { X -= C; }
                    }
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
                cblas_dsyrk(CblasColMajor,CblasUpper,CblasTrans,(int)C,(int)R,1.0,X,(int)R,0.0,Y,(int)C);
                for (int c=0; c<(int)C-1; ++c)
                {
                    Y += c + 1;
                    for (int r=c+1; r<(int)C; ++r, ++Y) { *Y = *(Y+(int)C*(r-c)+c-r); }
                }
            }
            else
            {
                cblas_dsyrk(CblasRowMajor,CblasLower,CblasTrans,(int)C,(int)R,1.0,X,(int)C,0.0,Y,(int)C);
                for (int r=0; r<(int)C-1; ++r)
                {
                    Y += r + 1;
                    for (int c=r+1; c<(int)C; ++c, ++Y) { *Y = *(Y+(int)C*(c-r)+r-c); }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                cblas_dsyrk(CblasColMajor,CblasUpper,CblasNoTrans,(int)R,(int)C,1.0,X,(int)R,0.0,Y,(int)R);
                for (int c=0; c<(int)R-1; ++c)
                {
                    Y += c + 1;
                    for (int r=c+1; r<(int)R; ++r, ++Y) { *Y = *(Y+(int)R*(r-c)+c-r); }
                }
            }
            else
            {
                cblas_dsyrk(CblasRowMajor,CblasLower,CblasNoTrans,(int)R,(int)C,1.0,X,(int)C,0.0,Y,(int)R);
                for (int r=0; r<(int)R-1; ++r)
                {
                    Y += r + 1;
                    for (int c=r+1; c<(int)R; ++c, ++Y) { *Y = *(Y+(int)R*(c-r)+r-c); }
                }
            }
        }
    }

    return 0;
}


int matmul1t_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const char tr)
{
    const size_t N = R*C;
    float x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0u) {}
    else
    {
        if (tr)
        {
            if (iscolmajor)
            {
                for (size_t c2=0u; c2<C; ++c2)
                {
                    Y += 2u*c2;
                    sm2r = 0.0f;
                    for (size_t r=0u; r<R; ++r, ++X)
                    {
                        x1r = x2r = *X; ++X;
                        x1i = -*X; x2i = *X;
                        sm2r += x1r*x2r - x1i*x2i;
                    }
                    *Y++ = sm2r; *Y++ = 0.0f;
                    X -= 2u*R;
                    for (size_t c1=c2+1; c1<C; ++c1, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t r=0u; r<R; ++r, ++X)
                        {
                            x1r = *X; x2r = *(X+2u*R*(c1-c2));
                            ++X;
                            x1i = -*X; x2i = *(X+2u*R*(c1-c2));
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = *(Y+2*((int)((c1-c2)*C+c2)-(int)c1)) = sm2r;
                        ++Y;
                        *Y = -sm2i; *(Y+2*((int)((c1-c2)*C+c2)-(int)c1)) = sm2i;
                        if (c1<C-1u) { X -= 2u*R; }
                    }
                }
            }
            else
            {
                for (size_t c1=0u; c1<C; ++c1, X+=2)
                {
                    Y += 2u*c1;
                    sm2r = 0.0f;
                    for (size_t r=0u; r<R; ++r, X+=2u*C-1u)
                    {
                        x1r = x2r = *X; ++X;
                        x1i = -*X; x2i = *X;
                        sm2r += x1r*x2r - x1i*x2i;
                    }
                    *Y++ = sm2r; *Y++ = 0.0f; X -= 2u*N;
                    for (size_t c2=c1+1; c2<C; ++c2, X-=2u*N, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t r=0u; r<R; ++r, X+=2u*C-1u)
                        {
                            x1r = *X; x2r = *(X+2u*(c2-c1));
                            ++X;
                            x1i = -*X; x2i = *(X+2u*(c2-c1));
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = *(Y+2*((int)((c2-c1)*C+c1)-(int)c2)) = sm2r;
                        ++Y;
                        *Y = sm2i; *(Y+2*((int)((c2-c1)*C+c1)-(int)c2)) = -sm2i;
                    }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (size_t r2=0u; r2<R; ++r2, X+=2)
                {
                    Y += 2u*r2;
                    sm2r = 0.0f;
                    for (size_t c=0u; c<C; ++c, X+=2u*R-1u)
                    {
                        x1r = x2r = *X; ++X;
                        x1i = *X; x2i = -*X;
                        sm2r += x1r*x2r - x1i*x2i;
                    }
                    *Y++ = sm2r; *Y++ = 0.0f; X -= 2u*N;
                    for (size_t r1=r2+1; r1<R; ++r1, X-=2u*N, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t c=0u; c<C; ++c, X+=2u*R-1u)
                        {
                            x1r = *X; x2r = *(X+2u*(r1-r2));
                            ++X;
                            x1i = *X; x2i = -*(X+2u*(r1-r2));
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = *(Y+2*((int)((r1-r2)*R+r2)-(int)r1)) = sm2r;
                        ++Y;
                        *Y = -sm2i; *(Y+2*((int)((r1-r2)*R+r2)-(int)r1)) = sm2i;
                    }
                }
            }
            else
            {
                for (size_t r1=0u; r1<R; ++r1)
                {
                    Y += 2u*r1;
                    sm2r = 0.0f;
                    for (size_t c=0u; c<C; ++c, ++X)
                    {
                        x1r = x2r = *X; ++X;
                        x1i = *X; x2i = -*X;
                        sm2r += x1r*x2r - x1i*x2i;
                    }
                    *Y++ = sm2r; *Y++ = 0.0f;
                    X -= 2u*C;
                    for (size_t r2=r1+1u; r2<R; ++r2, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t c=0u; c<C; ++c, ++X)
                        {
                            x1r = *X; x2r = *(X+2u*C*(r2-r1));
                            ++X;
                            x1i = *X; x2i = -*(X+2u*C*(r2-r1));
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = *(Y+2*((int)((r2-r1)*R+r1)-(int)r2)) = sm2r;
                        ++Y;
                        *Y = sm2i; *(Y+2*((int)((r2-r1)*R+r1)-(int)r2)) = -sm2i;
                        if (r2<R-1u) { X -= 2u*C; }
                    }
                }
            }
        }
    }

    return 0;
}


int matmul1t_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const char tr)
{
    const size_t N = R*C;
    double x1r, x1i, x2r, x2i, sm2r, sm2i;

    if (N==0u) {}
    else
    {
        if (tr)
        {
            if (iscolmajor)
            {
                for (size_t c2=0u; c2<C; ++c2)
                {
                    Y += 2u*c2;
                    sm2r = 0.0;
                    for (size_t r=0u; r<R; ++r, ++X)
                    {
                        x1r = x2r = *X; ++X;
                        x1i = -*X; x2i = *X;
                        sm2r += x1r*x2r - x1i*x2i;
                    }
                    *Y++ = sm2r; *Y++ = 0.0;
                    X -= 2u*R;
                    for (size_t c1=c2+1; c1<C; ++c1, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t r=0u; r<R; ++r, ++X)
                        {
                            x1r = *X; x2r = *(X+2u*R*(c1-c2));
                            ++X;
                            x1i = -*X; x2i = *(X+2u*R*(c1-c2));
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = *(Y+2*((int)((c1-c2)*C+c2)-(int)c1)) = sm2r;
                        ++Y;
                        *Y = -sm2i; *(Y+2*((int)((c1-c2)*C+c2)-(int)c1)) = sm2i;
                        if (c1<C-1u) { X -= 2u*R; }
                    }
                }
            }
            else
            {
                for (size_t c1=0u; c1<C; ++c1, X+=2)
                {
                    Y += 2u*c1;
                    sm2r = 0.0;
                    for (size_t r=0u; r<R; ++r, X+=2u*C-1u)
                    {
                        x1r = x2r = *X; ++X;
                        x1i = -*X; x2i = *X;
                        sm2r += x1r*x2r - x1i*x2i;
                    }
                    *Y++ = sm2r; *Y++ = 0.0; X -= 2u*N;
                    for (size_t c2=c1+1; c2<C; ++c2, X-=2u*N, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t r=0u; r<R; ++r, X+=2u*C-1u)
                        {
                            x1r = *X; x2r = *(X+2u*(c2-c1));
                            ++X;
                            x1i = -*X; x2i = *(X+2u*(c2-c1));
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = *(Y+2*((int)((c2-c1)*C+c1)-(int)c2)) = sm2r;
                        ++Y;
                        *Y = sm2i; *(Y+2*((int)((c2-c1)*C+c1)-(int)c2)) = -sm2i;
                    }
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (size_t r2=0u; r2<R; ++r2, X+=2)
                {
                    Y += 2u*r2;
                    sm2r = 0.0;
                    for (size_t c=0u; c<C; ++c, X+=2u*R-1u)
                    {
                        x1r = x2r = *X; ++X;
                        x1i = *X; x2i = -*X;
                        sm2r += x1r*x2r - x1i*x2i;
                    }
                    *Y++ = sm2r; *Y++ = 0.0; X -= 2u*N;
                    for (size_t r1=r2+1; r1<R; ++r1, X-=2u*N, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t c=0u; c<C; ++c, X+=2u*R-1u)
                        {
                            x1r = *X; x2r = *(X+2u*(r1-r2));
                            ++X;
                            x1i = *X; x2i = -*(X+2u*(r1-r2));
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = *(Y+2*((int)((r1-r2)*R+r2)-(int)r1)) = sm2r;
                        ++Y;
                        *Y = -sm2i; *(Y+2*((int)((r1-r2)*R+r2)-(int)r1)) = sm2i;
                    }
                }
            }
            else
            {
                for (size_t r1=0u; r1<R; ++r1)
                {
                    Y += 2u*r1;
                    sm2r = 0.0;
                    for (size_t c=0u; c<C; ++c, ++X)
                    {
                        x1r = x2r = *X; ++X;
                        x1i = *X; x2i = -*X;
                        sm2r += x1r*x2r - x1i*x2i;
                    }
                    *Y++ = sm2r; *Y++ = 0.0;
                    X -= 2u*C;
                    for (size_t r2=r1+1u; r2<R; ++r2, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t c=0u; c<C; ++c, ++X)
                        {
                            x1r = *X; x2r = *(X+2u*C*(r2-r1));
                            ++X;
                            x1i = *X; x2i = -*(X+2u*C*(r2-r1));
                            sm2r += x1r*x2r - x1i*x2i;
                            sm2i += x1r*x2i + x1i*x2r;
                        }
                        *Y = *(Y+2*((int)((r2-r1)*R+r1)-(int)r2)) = sm2r;
                        ++Y;
                        *Y = sm2i; *(Y+2*((int)((r2-r1)*R+r1)-(int)r2)) = -sm2i;
                        if (r2<R-1u) { X -= 2u*C; }
                    }
                }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
