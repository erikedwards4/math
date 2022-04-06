//Linear algebra function.
//Affine transformation of each vector in X using matrix A and vector B.

//This is different overall setup than usual matmul2: X can be 4D tensor,
//with the vecs along any dimension, and each is transformed (vec2vec operation).
//Thus, the code is much like the code for other vec2vec operations.

//Each vector in X has length Lx, and each vector in Y has length Ly.
//Vector B always has length Ly.
//If colmajor, then A has size Lx x Ly, and y = A'*x + B for each col vec x.
//If rowmajor, then A has size Ly x Lx, and y = A *x + B for each col vec x.
//Or:
//If colmajor, then A has size Lx x Ly, and y = x*A  + B for each row vec x.
//If rowmajor, then A has size Ly x Lx, and y = x*A' + B for each row vec x.
//That is:
//For performance reasons, this assumes that A has leading dimension Lx!

//Note: using cblas_?gemm is not faster (and slower at small N), becuase B must be copied into Y.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int affine_s (float *Y, const float *X, const float *A, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const int iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Na = Lx*Ly;
    float sm2;

    if (N==0u) {}
    else if (Lx==N)
    {
        if (Na<30000u)
        {
            for (size_t ly=Ly; ly>0u; --ly, X-=Lx, ++B, ++Y)
            {
                sm2 = *B;
                for (size_t lx=Lx; lx>0u; --lx, ++X, ++A) { sm2 = fmaf(*X,*A,sm2); }
                *Y = sm2;
            }
        }
        else
        {
            for (size_t ly=Ly; ly>0u; --ly, ++B, ++Y) { *Y = *B; }
            Y -= Ly;
            cblas_sgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,1.0f,A,(int)Lx,X,1,1.0f,Y,1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            if (Na<1000u)
            {
                for (size_t v=V; v>0u; --v, A-=Lx*Ly, B-=Ly)
                {
                    for (size_t ly=Ly; ly>0u; --ly, ++B, ++Y)
                    {
                        sm2 = *B;
                        for (size_t lx=Lx; lx>0u; --lx, ++X, ++A) { sm2 = fmaf(*X,*A,sm2); }
                        *Y = sm2;
                        if (ly>1u) { X -= Lx; }
                    }
                }
            }
            else
            {
                for (size_t v=V; v>0u; --v, X+=Lx, B-=Ly, Y+=Ly)
                {
                    for (size_t ly=Ly; ly>0u; --ly, ++B, ++Y) { *Y = *B; }
                    Y -= Ly;
                    cblas_sgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,1.0f,A,(int)Lx,X,1,1.0f,Y,1);
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=BS*(Lx-1u), Y+=BS*(Ly-1u))
            {
                for (size_t b=0u; b<BS; ++b, ++X, A-=Lx*Ly, B-=Ly, Y-=K*Ly-1u)
                {
                    for (size_t ly=Ly; ly>0u; --ly, X-=K*Lx, ++B, Y+=K)
                    {
                        sm2 = *B;
                        for (size_t lx=Lx; lx>0u; --lx, X+=K, ++A) { sm2 = fmaf(*X,*A,sm2); }
                        *Y = sm2;
                    }
                }
            }
        }
    }

    return 0;
}


int affine_d (double *Y, const double *X, const double *A, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const int iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Na = Lx*Ly;
    double sm2;

    if (N==0u) {}
    else if (Lx==N)
    {
        if (Na<30000u)
        {
            for (size_t ly=Ly; ly>0u; --ly, X-=Lx, ++B, ++Y)
            {
                sm2 = *B;
                for (size_t lx=Lx; lx>0u; --lx, ++X, ++A) { sm2 = fma(*X,*A,sm2); }
                *Y = sm2;
            }
        }
        else
        {
            for (size_t ly=Ly; ly>0u; --ly, ++B, ++Y) { *Y = *B; }
            Y -= Ly;
            cblas_dgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,1.0,A,(int)Lx,X,1,1.0,Y,1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            if (Na<1000u)
            {
                for (size_t v=V; v>0u; --v, A-=Lx*Ly, B-=Ly)
                {
                    for (size_t ly=Ly; ly>0u; --ly, ++B, ++Y)
                    {
                        sm2 = *B;
                        for (size_t lx=Lx; lx>0u; --lx, ++X, ++A) { sm2 = fma(*X,*A,sm2); }
                        *Y = sm2;
                        if (ly>1u) { X -= Lx; }
                    }
                }
            }
            else
            {
                for (size_t v=V; v>0u; --v, X+=Lx, B-=Ly, Y+=Ly)
                {
                    for (size_t ly=Ly; ly>0u; --ly, ++B, ++Y) { *Y = *B; }
                    Y -= Ly;
                    cblas_dgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,1.0,A,(int)Lx,X,1,1.0,Y,1);
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=BS*(Lx-1u), Y+=BS*(Ly-1u))
            {
                for (size_t b=0u; b<BS; ++b, ++X, A-=Lx*Ly, B-=Ly, Y-=K*Ly-1u)
                {
                    for (size_t ly=Ly; ly>0u; --ly, X-=K*Lx, ++B, Y+=K)
                    {
                        sm2 = *B;
                        for (size_t lx=Lx; lx>0u; --lx, X+=K, ++A) { sm2 = fma(*X,*A,sm2); }
                        *Y = sm2;
                    }
                }
            }
        }
    }

    return 0;
}


int affine_c (float *Y, const float *X, const float *A, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const int iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Na = Lx*Ly;
    float xr, xi, ar, ai, sm2r, sm2i;

    if (N==0u) {}
    else if (Lx==N)
    {
        if (Lx<30000)
        {
            for (size_t ly=Ly; ly>0u; --ly, X-=2u*Lx, ++B, ++Y)
            {
                sm2r = *B; sm2i = *++B;
                for (size_t lx=Lx; lx>0u; --lx, ++X, ++A)
                {
                    xr = *X; xi = *++X;
                    ar = *A; ai = *++A;
                    sm2r += xr*ar - xi*ai;
                    sm2i += xr*ai + xi*ar;
                }
                *Y = sm2r; *++Y = sm2i;
            }
        }
        else
        {
            const float o[2] = {1.0f,0.0f};
            for (size_t ly=2u*Ly; ly>0u; --ly, ++B, ++Y) { *Y = *B; }
            Y -= 2u*Ly;
            cblas_cgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,o,A,(int)Lx,X,1,o,Y,1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            if (Na<4500u)
            {
                for (size_t v=V; v>0u; --v, A-=2u*Na, B-=2u*Ly)
                {
                    for (size_t ly=Ly; ly>0u; --ly, ++B, ++Y)
                    {
                        sm2r = *B; sm2i = *++B;
                        for (size_t lx=Lx; lx>0u; --lx, ++X, ++A)
                        {
                            xr = *X; xi = *++X;
                            ar = *A; ai = *++A;
                            sm2r += xr*ar - xi*ai;
                            sm2i += xr*ai + xi*ar;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (ly>1u) { X -= 2u*Lx; }
                    }
                }
            }
            else
            {
                const float o[2] = {1.0f,0.0f};
                for (size_t v=V; v>0u; --v, X+=2u*Lx, B-=2u*Ly, Y+=2u*Ly)
                {
                    for (size_t ly=2u*Ly; ly>0u; --ly, ++B, ++Y) { *Y = *B; }
                    Y -= 2u*Ly;
                    cblas_cgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,o,A,(int)Lx,X,1,o,Y,1);
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*BS*(Lx-1u), Y+=2u*BS*(Ly-1u))
            {
                for (size_t b=0u; b<BS; ++b, X+=2, A-=2u*Na, B-=2u*Ly, Y-=2u*K*Ly-2u)
                {
                    for (size_t ly=Ly; ly>0u; --ly, X-=2u*K*Lx, ++B, Y+=2u*K-1u)
                    {
                        sm2r = *B; sm2i = *++B;
                        for (size_t lx=Lx; lx>0u; --lx, X+=2u*K-1u, ++A)
                        {
                            xr = *X; xi = *++X;
                            ar = *A; ai = *++A;
                            sm2r += xr*ar - xi*ai;
                            sm2i += xr*ai + xi*ar;
                        }
                        *Y = sm2r; *++Y = sm2i;
                    }
                }
            }
        }
    }

    return 0;
}


int affine_z (double *Y, const double *X, const double *A, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const int iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Na = Lx*Ly;
    double xr, xi, ar, ai, sm2r, sm2i;

    if (N==0u) {}
    else if (Lx==N)
    {
        if (Lx<30000)
        {
            for (size_t ly=Ly; ly>0u; --ly, X-=2u*Lx, ++B, ++Y)
            {
                sm2r = *B; sm2i = *++B;
                for (size_t lx=Lx; lx>0u; --lx, ++X, ++A)
                {
                    xr = *X; xi = *++X;
                    ar = *A; ai = *++A;
                    sm2r += xr*ar - xi*ai;
                    sm2i += xr*ai + xi*ar;
                }
                *Y = sm2r; *++Y = sm2i;
            }
        }
        else
        {
            const double o[2] = {1.0,0.0};
            for (size_t ly=2u*Ly; ly>0u; --ly, ++B, ++Y) { *Y = *B; }
            Y -= 2u*Ly;
            cblas_zgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,o,A,(int)Lx,X,1,o,Y,1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            if (Na<4500u)
            {
                for (size_t v=V; v>0u; --v, A-=2u*Na, B-=2u*Ly)
                {
                    for (size_t ly=Ly; ly>0u; --ly, ++B, ++Y)
                    {
                        sm2r = *B; sm2i = *++B;
                        for (size_t lx=Lx; lx>0u; --lx, ++X, ++A)
                        {
                            xr = *X; xi = *++X;
                            ar = *A; ai = *++A;
                            sm2r += xr*ar - xi*ai;
                            sm2i += xr*ai + xi*ar;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (ly>1u) { X -= 2u*Lx; }
                    }
                }
            }
            else
            {
                const double o[2] = {1.0,0.0};
                for (size_t v=V; v>0u; --v, X+=2u*Lx, B-=2u*Ly, Y+=2u*Ly)
                {
                    for (size_t ly=2u*Ly; ly>0u; --ly, ++B, ++Y) { *Y = *B; }
                    Y -= 2u*Ly;
                    cblas_zgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,o,A,(int)Lx,X,1,o,Y,1);
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*BS*(Lx-1u), Y+=2u*BS*(Ly-1u))
            {
                for (size_t b=0u; b<BS; ++b, X+=2, A-=2u*Na, B-=2u*Ly, Y-=2u*K*Ly-2u)
                {
                    for (size_t ly=Ly; ly>0u; --ly, X-=2u*K*Lx, ++B, Y+=2u*K-1u)
                    {
                        sm2r = *B; sm2i = *++B;
                        for (size_t lx=Lx; lx>0u; --lx, X+=2u*K-1u, ++A)
                        {
                            xr = *X; xi = *++X;
                            ar = *A; ai = *++A;
                            sm2r += xr*ar - xi*ai;
                            sm2i += xr*ai + xi*ar;
                        }
                        *Y = sm2r; *++Y = sm2i;
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
