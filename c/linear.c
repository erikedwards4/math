//Linear algebra function.
//Linear transformation of each vector in X using matrix A.

//This is different overall setup than usual matmul2: X can be 4D tensor,
//with the vecs along any dimension, and each is transformed (vec2vec operation).
//Thus, the code is much like the code for other vec2vec operations.

//Each vector in X has length Lx, and each vector in Y has length Ly.
//If colmajor, then A has size Lx x Ly, and y = A'*x for each col vec x.
//If rowmajor, then A has size Ly x Lx, and y = A*x for each col vec x.
//Or:
//If colmajor, then A has size Lx x Ly, and y = x*A for each row vec x.
//If rowmajor, then A has size Ly x Lx, and y = x*A' for each row vec x.
//That is:
//For performance reasons, this assumes that A has leading dimension Lx!

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int linear_s (float *Y, const float *X, const float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim);
int linear_d (double *Y, const double *X, const double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim);
int linear_c (float *Y, const float *X, const float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim);
int linear_z (double *Y, const double *X, const double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim);


int linear_s (float *Y, const float *X, const float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim)
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
            for (size_t ly=0u; ly<Ly; ++ly, X-=Lx, ++Y)
            {
                sm2 = 0.0f;
                for (size_t lx=0u; lx<Lx; ++lx, ++X, ++A) { sm2 = fmaf(*X,*A,sm2); }
                *Y = sm2;
            }
        }
        else
        {
            cblas_sgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,1.0f,A,(int)Lx,X,1,0.0f,Y,1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (Na<4500u)
            {
                for (size_t v=0u; v<V; ++v, A-=Na)
                {
                    for (size_t ly=0u; ly<Ly; ++ly, ++Y)
                    {
                        sm2 = 0.0f;
                        for (size_t lx=0u; lx<Lx; ++lx, ++X, ++A) { sm2 = fmaf(*X,*A,sm2); }
                        *Y = sm2;
                        if (ly<Ly-1u) { X -= Lx; }
                    }
                }
            }
            else
            {
                cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,(int)Ly,(int)V,(int)Lx,1.0f,A,(int)Lx,X,(int)Lx,0.0f,Y,(int)Ly);
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=0u; b<B; ++b, ++X, A-=Na, Y-=K*Ly-1u)
                {
                    for (size_t ly=0u; ly<Ly; ++ly, X-=K*Lx, Y+=K)
                    {
                        sm2 = 0.0f;
                        for (size_t lx=0u; lx<Lx; ++lx, X+=K, ++A) { sm2 = fmaf(*X,*A,sm2); }
                        *Y = sm2;
                    }
                }
            }
        }
    }

    return 0;
}


int linear_d (double *Y, const double *X, const double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim)
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
            for (size_t ly=0u; ly<Ly; ++ly, X-=Lx, ++Y)
            {
                sm2 = 0.0;
                for (size_t lx=0u; lx<Lx; ++lx, ++X, ++A) { sm2 = fma(*X,*A,sm2); }
                *Y = sm2;
            }
        }
        else
        {
            cblas_dgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,1.0,A,(int)Lx,X,1,0.0,Y,1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (Na<4500u)
            {
                for (size_t v=0u; v<V; ++v, A-=Na)
                {
                    for (size_t ly=0u; ly<Ly; ++ly, ++Y)
                    {
                        sm2 = 0.0;
                        for (size_t lx=0u; lx<Lx; ++lx, ++X, ++A) { sm2 = fma(*X,*A,sm2); }
                        *Y = sm2;
                        if (ly<Ly-1u) { X -= Lx; }
                    }
                }
            }
            else
            {
                cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,(int)Ly,(int)V,(int)Lx,1.0,A,(int)Lx,X,(int)Lx,0.0,Y,(int)Ly);
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=0u; b<B; ++b, ++X, A-=Na, Y-=K*Ly-1u)
                {
                    for (size_t ly=0u; ly<Ly; ++ly, X-=K*Lx, Y+=K)
                    {
                        sm2 = 0.0;
                        for (size_t lx=0u; lx<Lx; ++lx, X+=K, ++A) { sm2 = fma(*X,*A,sm2); }
                        *Y = sm2;
                    }
                }
            }
        }
    }

    return 0;
}


int linear_c (float *Y, const float *X, const float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Na = Lx*Ly;
    float sm2r, sm2i, xr, xi, ar, ai;

    if (N==0u) {}
    else if (Lx==N)
    {
        if (Na<30000u)
        {
            for (size_t ly=0u; ly<Ly; ++ly, X-=2u*Lx, ++Y)
            {
                sm2r = sm2i = 0.0f;
                for (size_t lx=0u; lx<Lx; ++lx, ++X, ++A)
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
            const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
            cblas_cgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,o,A,(int)Lx,X,1,z,Y,1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (Na<4500u)
            {
                for (size_t v=0u; v<V; ++v, A-=2u*Na)
                {
                    for (size_t ly=0u; ly<Ly; ++ly, ++Y)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t lx=0u; lx<Lx; ++lx, ++X, ++A)
                        {
                            xr = *X; xi = *++X;
                            ar = *A; ai = *++A;
                            sm2r += xr*ar - xi*ai;
                            sm2i += xr*ai + xi*ar;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (ly<Ly-1u) { X -= 2u*Lx; }
                    }
                }
            }
            else
            {
                const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
                cblas_cgemm(CblasColMajor,CblasTrans,CblasNoTrans,(int)Ly,(int)V,(int)Lx,o,A,(int)Lx,X,(int)Lx,z,Y,(int)Ly);
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
            {
                for (size_t b=0u; b<B; ++b, ++X, A-=2u*Na, Y-=2u*K*Ly-2u)
                {
                    for (size_t ly=0u; ly<Ly; ++ly, X-=2u*K*Lx, Y+=2u*K-1u)
                    {
                        sm2r = sm2i = 0.0f;
                        for (size_t lx=0u; lx<Lx; ++lx, X+=2u*K-1u, ++A)
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


int linear_z (double *Y, const double *X, const double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Ly, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Na = Lx*Ly;
    double sm2r, sm2i, xr, xi, ar, ai;

    if (N==0u) {}
    else if (Lx==N)
    {
        if (Na<150)
        {
            for (size_t ly=0u; ly<Ly; ++ly, X-=2u*Lx, ++Y)
            {
                sm2r = sm2i = 0.0;
                for (size_t lx=0u; lx<Lx; ++lx, ++X, ++A)
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
            const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
            cblas_zgemv(CblasColMajor,CblasTrans,(int)Lx,(int)Ly,o,A,(int)Lx,X,1,z,Y,1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (Na<4500u)
            {
                for (size_t v=0u; v<V; ++v, A-=2u*Na)
                {
                    for (size_t ly=0u; ly<Ly; ++ly, ++Y)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t lx=0u; lx<Lx; ++lx, ++X, ++A)
                        {
                            xr = *X; xi = *++X;
                            ar = *A; ai = *++A;
                            sm2r += xr*ar - xi*ai;
                            sm2i += xr*ai + xi*ar;
                        }
                        *Y = sm2r; *++Y = sm2i;
                        if (ly<Ly-1u) { X -= 2u*Lx; }
                    }
                }
            }
            else
            {
                const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
                cblas_zgemm(CblasColMajor,CblasTrans,CblasNoTrans,(int)Ly,(int)V,(int)Lx,o,A,(int)Lx,X,(int)Lx,z,Y,(int)Ly);
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
            {
                for (size_t b=0u; b<B; ++b, ++X, A-=2u*Na, Y-=2u*K*Ly-2u)
                {
                    for (size_t ly=0u; ly<Ly; ++ly, X-=2u*K*Lx, Y+=2u*K-1u)
                    {
                        sm2r = sm2i = 0.0;
                        for (size_t lx=0u; lx<Lx; ++lx, X+=2u*K-1u, ++A)
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
