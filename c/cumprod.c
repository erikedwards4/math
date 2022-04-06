//Cumulative product along each vector in X (vec2vec operation)
//This has in-place and not-in-place versions.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int cumprod_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cumprod_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y++ = *X++;
        for (size_t l=L; l>1u; --l, ++X, ++Y) { *Y = *(Y-1) * *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                *Y++ = *X++;
                for (size_t l=L; l>1u; --l, ++X, ++Y) { *Y = *(Y-1) * *X; }
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; }
            for (size_t l=L; l>1u; --l)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *(Y-V) * *X; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                {
                    *Y = *X; X += K; Y += K;
                    for (size_t l=L; l>1u; --l, X+=K, Y+=K) { *Y = *(Y-K) * *X; }
                }
            }
        }
    }

    return 0;
}


int cumprod_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cumprod_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y++ = *X++;
        for (size_t l=L; l>1u; --l, ++X, ++Y) { *Y = *(Y-1) * *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                *Y++ = *X++;
                for (size_t l=L; l>1u; --l, ++X, ++Y) { *Y = *(Y-1) * *X; }
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; }
            for (size_t l=L; l>1u; --l)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *(Y-V) * *X; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                {
                    *Y = *X; X += K; Y += K;
                    for (size_t l=L; l>1u; --l, X+=K, Y+=K) { *Y = *(Y-K) * *X; }
                }
            }
        }
    }

    return 0;
}


int cumprod_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cumprod_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y++ = *X++; *Y++ = *X++;
        for (size_t l=L; l>1u; --l, X+=2, Y+=2)
        {
            *Y = *(Y-2)**X - *(Y-1)**(X+1);
            *(Y+1) = *(Y-1)**X + *(Y-2)**(X+1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                *Y++ = *X++; *Y++ = *X++;
                for (size_t l=L; l>1u; --l, ++X, ++Y) { *Y = *(Y-2) + *X; ++Y; *Y = *(Y-2) + *++X; }
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            for (size_t l=L; l>1u; --l)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *(Y-2u*V) + *X; ++Y; *Y = *(Y-2u*V) + *++X; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y-=2u*K*L-2u)
                {
                    *Y = *X; *(Y+1) = *(X+1); X += 2u*K; Y += 2u*K;
                    for (size_t l=L; l>1u; --l, X+=2u*K, Y+=2u*K) { *Y = *(Y-2u*K) + *X; *(Y+1) = *(Y-2u*K+1u) + *(X+1); }
                }
            }
        }
    }

    return 0;
}


int cumprod_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cumprod_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y++ = *X++; *Y++ = *X++;
        for (size_t l=L; l>1u; --l, ++Y) { *Y = *(Y-2) + *X++; ++Y; *Y = *(Y-2) + *X++; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                *Y++ = *X++; *Y++ = *X++;
                for (size_t l=L; l>1u; --l, ++X, ++Y) { *Y = *(Y-2) + *X; ++Y; *Y = *(Y-2) + *++X; }
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            for (size_t l=L; l>1u; --l)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *(Y-2u*V) + *X; ++Y; *Y = *(Y-2u*V) + *++X; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y-=2u*K*L-2u)
                {
                    *Y = *X; *(Y+1) = *(X+1); X += 2u*K; Y += 2u*K;
                    for (size_t l=L; l>1u; --l, X+=2u*K, Y+=2u*K) { *Y = *(Y-2u*K) + *X; *(Y+1) = *(Y-2u*K+1u) + *(X+1); }
                }
            }
        }
    }

    return 0;
}


int cumprod_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cumprod_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        ++X;
        for (size_t l=L; l>1u; --l, ++X) { *X *= *(X-1); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                ++X;
                for (size_t l=L; l>1u; --l, ++X) { *X *= *(X-1); }
            }
        }
        else if (G==1u)
        {
            X += V;
            for (size_t l=L; l>1u; --l)
            {
                for (size_t v=V; v>0u; --v, ++X) { *X *= *(X-V); }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u)
                {
                    X += K;
                    for (size_t l=L; l>1u; --l, X+=K) { *X *= *(X-K); }
                }
            }
        }
    }

    return 0;
}


int cumprod_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cumprod_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        ++X;
        for (size_t l=L; l>1u; --l, ++X) { *X *= *(X-1); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                ++X;
                for (size_t l=L; l>1u; --l, ++X) { *X *= *(X-1); }
            }
        }
        else if (G==1u)
        {
            X += V;
            for (size_t l=L; l>1u; --l)
            {
                for (size_t v=V; v>0u; --v, ++X) { *X *= *(X-V); }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u)
                {
                    X += K;
                    for (size_t l=L; l>1u; --l, X+=K) { *X *= *(X-K); }
                }
            }
        }
    }

    return 0;
}


int cumprod_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cumprod_inplace_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float xr, xi, c, d;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        c = *X++; d = *X++;
        for (size_t l=L; l>1u; --l, ++X)
        {
            xr = *X; xi = *(X+1);
            *X = xr*c - xi*d;
            *++X = xr*d + xi*c;
            c = *(X-1); d = *X;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                c = *X++; d = *X++;
                for (size_t l=L; l>1u; --l, ++X)
                {
                    xr = *X; xi = *(X+1);
                    *X = xr*c - xi*d;
                    *++X = xr*d + xi*c;
                    c = *(X-1); d = *X;
                }
            }
        }
        else if (G==1u)
        {
            X += 2u*V;
            for (size_t l=L; l>1u; --l)
            {
                for (size_t v=V; v>0u; --v, ++X)
                {
                    c = *(X-2u*V); d = *(X-2u*V+1);
                    xr = *X; xi = *(X+1);
                    *X = xr*c - xi*d;
                    *++X = xr*d + xi*c;
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u)
                {
                    c = *X; d = *++X; X += 2u*K-1u;
                    for (size_t l=L; l>1u; --l, X+=2u*K-1u)
                    {
                        xr = *X; xi = *(X+1);
                        *X = xr*c - xi*d;
                        *++X = xr*d + xi*c;
                        c = *(X-1); d = *X;
                    }
                }
            }
        }
    }

    return 0;
}


int cumprod_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cumprod_inplace_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double xr, xi, c, d;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        c = *X++; d = *X++;
        for (size_t l=L; l>1u; --l, ++X)
        {
            xr = *X; xi = *(X+1);
            *X = xr*c - xi*d;
            *++X = xr*d + xi*c;
            c = *(X-1); d = *X;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                c = *X++; d = *X++;
                for (size_t l=L; l>1u; --l, ++X)
                {
                    xr = *X; xi = *(X+1);
                    *X = xr*c - xi*d;
                    *++X = xr*d + xi*c;
                    c = *(X-1); d = *X;
                }
            }
        }
        else if (G==1u)
        {
            X += 2u*V;
            for (size_t l=L; l>1u; --l)
            {
                for (size_t v=V; v>0u; --v, ++X)
                {
                    c = *(X-2u*V); d = *(X-2u*V+1);
                    xr = *X; xi = *(X+1);
                    *X = xr*c - xi*d;
                    *++X = xr*d + xi*c;
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u)
                {
                    c = *X; d = *++X; X += 2u*K-1u;
                    for (size_t l=L; l>1u; --l, X+=2u*K-1u)
                    {
                        xr = *X; xi = *(X+1);
                        *X = xr*c - xi*d;
                        *++X = xr*d + xi*c;
                        c = *(X-1); d = *X;
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
