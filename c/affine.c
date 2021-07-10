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

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int affine_s (float *Y, const float *X, const float *A, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim);
int affine_d (double *Y, const double *X, const double *A, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim);
int affine_c (float *Y, const float *X, const float *A, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim);
int affine_z (double *Y, const double *X, const double *A, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim);


int affine_s (float *Y, const float *X, const float *A, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Na = Lx*Ly;
    float sm2;

    if (N==0u) {}
    else if (Lx==N)
    {
        if (Na<30000)
        {
            for (size_t ly=0; ly<Ly; ++ly, X-=Lx, ++B, ++Y)
            {
                sm2 = *B;
                for (size_t lx=0; lx<Lx; ++lx, ++X, ++A) { sm2 = fmaf(*X,*A,sm2); }
                *Y = sm2;
            }
        }
        else
        {
            for (size_t ly=0; ly<Ly; ++ly, ++B, ++Y) { *Y = *B; }
            Y -= Ly;
            cblas_sgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,1.0f,A,(int)Lx,X,1,1.0f,Y,1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t BS = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;

        if (K==1 && (G==1 || BS==1))
        {
            if (Na<1000)
            {
                for (size_t v=0u; v<V; ++v, A-=Lx*Ly, B-=Ly)
                {
                    for (size_t ly=0; ly<Ly; ++ly, ++B, ++Y)
                    {
                        sm2 = *B;
                        for (size_t lx=0; lx<Lx; ++lx, ++X, ++A) { sm2 = fmaf(*X,*A,sm2); }
                        *Y = sm2;
                        if (ly<Ly-1) { X -= Lx; }
                    }
                }
            }
            else
            {
                for (size_t v=0u; v<V; ++v, X+=Lx, B-=Ly, Y+=Ly)
                {
                    for (size_t ly=0; ly<Ly; ++ly, ++B, ++Y) { *Y = *B; }
                    Y -= Ly;
                    cblas_sgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,1.0f,A,(int)Lx,X,1,1.0f,Y,1);
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=BS*(Lx-1), Y+=BS*(Ly-1))
            {
                for (size_t b=0u; b<BS; ++b, ++X, A-=Lx*Ly, B-=Ly, Y-=K*Ly-1)
                {
                    for (size_t ly=0; ly<Ly; ++ly, X-=K*Lx, ++B, Y+=K)
                    {
                        sm2 = *B;
                        for (size_t lx=0; lx<Lx; ++lx, X+=K, ++A) { sm2 = fmaf(*X,*A,sm2); }
                        *Y = sm2;
                    }
                }
            }
        }
    }

    return 0;
}


int affine_d (double *Y, const double *X, const double *A, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Na = Lx*Ly;
    double sm2;

    if (N==0u) {}
    else if (Lx==N)
    {
        if (Na<30000)
        {
            for (size_t ly=0; ly<Ly; ++ly, X-=Lx, ++B, ++Y)
            {
                sm2 = *B;
                for (size_t lx=0; lx<Lx; ++lx, ++X, ++A) { sm2 = fma(*X,*A,sm2); }
                *Y = sm2;
            }
        }
        else
        {
            for (size_t ly=0; ly<Ly; ++ly, ++B, ++Y) { *Y = *B; }
            Y -= Ly;
            cblas_dgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,1.0,A,(int)Lx,X,1,1.0,Y,1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t BS = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;

        if (K==1 && (G==1 || BS==1))
        {
            if (Na<1000)
            {
                for (size_t v=0u; v<V; ++v, A-=Lx*Ly, B-=Ly)
                {
                    for (size_t ly=0; ly<Ly; ++ly, ++B, ++Y)
                    {
                        sm2 = *B;
                        for (size_t lx=0; lx<Lx; ++lx, ++X, ++A) { sm2 = fma(*X,*A,sm2); }
                        *Y = sm2;
                        if (ly<Ly-1) { X -= Lx; }
                    }
                }
            }
            else
            {
                for (size_t v=0u; v<V; ++v, X+=Lx, B-=Ly, Y+=Ly)
                {
                    for (size_t ly=0; ly<Ly; ++ly, ++B, ++Y) { *Y = *B; }
                    Y -= Ly;
                    cblas_dgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,1.0,A,(int)Lx,X,1,1.0,Y,1);
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=BS*(Lx-1), Y+=BS*(Ly-1))
            {
                for (size_t b=0u; b<BS; ++b, ++X, A-=Lx*Ly, B-=Ly, Y-=K*Ly-1)
                {
                    for (size_t ly=0; ly<Ly; ++ly, X-=K*Lx, ++B, Y+=K)
                    {
                        sm2 = *B;
                        for (size_t lx=0; lx<Lx; ++lx, X+=K, ++A) { sm2 = fma(*X,*A,sm2); }
                        *Y = sm2;
                    }
                }
            }
        }
    }

    return 0;
}


int affine_c (float *Y, const float *X, const float *A, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Na = Lx*Ly;
    float xr, xi, ar, ai, sm2r, sm2i;

    if (N==0u) {}
    else if (Lx==N)
    {
        if (Lx<30000)
        {
            for (size_t ly=0; ly<Ly; ++ly, X-=2*Lx, ++B, ++Y)
            {
                sm2r = *B; sm2i = *++B;
                for (size_t lx=0; lx<Lx; ++lx, ++X, ++A)
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
            for (size_t ly=0; ly<2*Ly; ++ly, ++B, ++Y) { *Y = *B; }
            Y -= 2*Ly;
            cblas_cgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,o,A,(int)Lx,X,1,o,Y,1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t BS = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;

        if (K==1 && (G==1 || BS==1))
        {
            if (Na<4500)
            {
                for (size_t v=0u; v<V; ++v, A-=2*Na, B-=2*Ly)
                {
                    for (size_t ly=0; ly<Ly; ++ly, ++B, ++Y)
                    {
                        sm2r = *B; sm2i = *++B;
                        for (size_t lx=0; lx<Lx; ++lx, ++X, ++A)
                        {
                            xr = *X; xi = *++X;
                            ar = *A; ai = *++A;
                            sm2r += xr*ar - xi*ai;
                            sm2i += xr*ai + xi*ar;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (ly<Ly-1) { X -= 2*Lx; }
                    }
                }
            }
            else
            {
                const float o[2] = {1.0f,0.0f};
                for (size_t v=0u; v<V; ++v, X+=2*Lx, B-=2*Ly, Y+=2*Ly)
                {
                    for (size_t ly=0; ly<2*Ly; ++ly, ++B, ++Y) { *Y = *B; }
                    Y -= 2*Ly;
                    cblas_cgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,o,A,(int)Lx,X,1,o,Y,1);
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*BS*(Lx-1), Y+=2*BS*(Ly-1))
            {
                for (size_t b=0u; b<BS; ++b, X+=2, A-=2*Na, B-=2*Ly, Y-=2*K*Ly-2)
                {
                    for (size_t ly=0; ly<Ly; ++ly, X-=2*K*Lx, ++B, Y+=2*K-1)
                    {
                        sm2r = *B; sm2i = *++B;
                        for (size_t lx=0; lx<Lx; ++lx, X+=2*K-1, ++A)
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


int affine_z (double *Y, const double *X, const double *A, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Na = Lx*Ly;
    double xr, xi, ar, ai, sm2r, sm2i;

    if (N==0u) {}
    else if (Lx==N)
    {
        if (Lx<30000)
        {
            for (size_t ly=0; ly<Ly; ++ly, X-=2*Lx, ++B, ++Y)
            {
                sm2r = *B; sm2i = *++B;
                for (size_t lx=0; lx<Lx; ++lx, ++X, ++A)
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
            for (size_t ly=0; ly<2*Ly; ++ly, ++B, ++Y) { *Y = *B; }
            Y -= 2*Ly;
            cblas_zgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,o,A,(int)Lx,X,1,o,Y,1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t BS = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/BS;

        if (K==1 && (G==1 || BS==1))
        {
            if (Na<4500)
            {
                for (size_t v=0u; v<V; ++v, A-=2*Na, B-=2*Ly)
                {
                    for (size_t ly=0; ly<Ly; ++ly, ++B, ++Y)
                    {
                        sm2r = *B; sm2i = *++B;
                        for (size_t lx=0; lx<Lx; ++lx, ++X, ++A)
                        {
                            xr = *X; xi = *++X;
                            ar = *A; ai = *++A;
                            sm2r += xr*ar - xi*ai;
                            sm2i += xr*ai + xi*ar;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (ly<Ly-1) { X -= 2*Lx; }
                    }
                }
            }
            else
            {
                const double o[2] = {1.0,0.0};
                for (size_t v=0u; v<V; ++v, X+=2*Lx, B-=2*Ly, Y+=2*Ly)
                {
                    for (size_t ly=0; ly<2*Ly; ++ly, ++B, ++Y) { *Y = *B; }
                    Y -= 2*Ly;
                    cblas_zgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,o,A,(int)Lx,X,1,o,Y,1);
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*BS*(Lx-1), Y+=2*BS*(Ly-1))
            {
                for (size_t b=0u; b<BS; ++b, X+=2, A-=2*Na, B-=2*Ly, Y-=2*K*Ly-2)
                {
                    for (size_t ly=0; ly<Ly; ++ly, X-=2*K*Lx, ++B, Y+=2*K-1)
                    {
                        sm2r = *B; sm2i = *++B;
                        for (size_t lx=0; lx<Lx; ++lx, X+=2*K-1, ++A)
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
