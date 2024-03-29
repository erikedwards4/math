//Vec2vec operation.
//Postpads each vector in X with P elements equal to val.
//Thus, Y has the same size as X along all other dims, but has a greater length than X along dim.

//For complex case, val is used for both real and imag parts.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int postpad_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const float val)
{
    if (dim>3u) { fprintf(stderr,"error in postpad_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u && Lx>0) {}
    else if (Lx==N)
    {
        for (size_t l=Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
        for (size_t p=P; p>0u; --p, ++Y) { *Y = val; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
                for (size_t p=P; p>0u; --p, ++Y) { *Y = val; }
            }
        }
        else if (G==1u)
        {
            for (size_t n=V*Lx; n>0u; --n, ++X, ++Y) { *Y = *X; }
            for (size_t n=V*P; n>0u; --n, ++Y) { *Y = val; }
        }
        else
        {
            const size_t Ly = Lx + P;
            for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=K*Ly-1u)
                {
                    for (size_t l=Lx; l>0u; --l, X+=K, Y+=K) { *Y = *X; }
                    for (size_t p=P; p>0u; --p, Y+=K) { *Y = val; }
                }
            }
        }
    }

    return 0;
}


int postpad_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const double val)
{
    if (dim>3u) { fprintf(stderr,"error in postpad_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u && Lx>0) {}
    else if (Lx==N)
    {
        for (size_t l=Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
        for (size_t p=P; p>0u; --p, ++Y) { *Y = val; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
                for (size_t p=P; p>0u; --p, ++Y) { *Y = val; }
            }
        }
        else if (G==1u)
        {
            for (size_t n=V*Lx; n>0u; --n, ++X, ++Y) { *Y = *X; }
            for (size_t n=V*P; n>0u; --n, ++Y) { *Y = val; }
        }
        else
        {
            const size_t Ly = Lx + P;
            for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=K*Ly-1u)
                {
                    for (size_t l=Lx; l>0u; --l, X+=K, Y+=K) { *Y = *X; }
                    for (size_t p=P; p>0u; --p, Y+=K) { *Y = val; }
                }
            }
        }
    }

    return 0;
}


int postpad_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const float val)
{
    if (dim>3u) { fprintf(stderr,"error in postpad_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u && Lx>0) {}
    else if (Lx==N)
    {
        for (size_t l=2u*Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
        for (size_t p=2u*P; p>0u; --p, ++Y) { *Y = val; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=2u*Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
                for (size_t p=2u*P; p>0u; --p, ++Y) { *Y = val; }
            }
        }
        else if (G==1u)
        {
            for (size_t n=2u*V*Lx; n>0u; --n, ++X, ++Y) { *Y = *X; }
            for (size_t n=2u*V*P; n>0u; --n, ++Y) { *Y = val; }
        }
        else
        {
            const size_t Ly = Lx + P;
            for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=2u*K*Ly-2u)
                {
                    for (size_t l=Lx; l>0u; --l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
                    for (size_t p=P; p>0u; --p, Y+=2u*K-1u) { *Y = val; *++Y = val; }
                }
            }
        }
    }

    return 0;
}


int postpad_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const double val)
{
    if (dim>3u) { fprintf(stderr,"error in postpad_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u && Lx>0) {}
    else if (Lx==N)
    {
        for (size_t l=2u*Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
        for (size_t p=2u*P; p>0u; --p, ++Y) { *Y = val; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=2u*Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
                for (size_t p=2u*P; p>0u; --p, ++Y) { *Y = val; }
            }
        }
        else if (G==1u)
        {
            for (size_t n=2u*V*Lx; n>0u; --n, ++X, ++Y) { *Y = *X; }
            for (size_t n=2u*V*P; n>0u; --n, ++Y) { *Y = val; }
        }
        else
        {
            const size_t Ly = Lx + P;
            for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=2u*K*Ly-2u)
                {
                    for (size_t l=Lx; l>0u; --l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
                    for (size_t p=P; p>0u; --p, Y+=2u*K-1u) { *Y = val; *++Y = val; }
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
