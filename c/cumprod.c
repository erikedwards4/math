//Cumulative product along each vector in X (vec2vec operation)
//This has in-place and not-in-place versions.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cumprod_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cumprod_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cumprod_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cumprod_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int cumprod_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cumprod_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cumprod_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cumprod_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int cumprod_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cumprod_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y++ = *X++;
        for (size_t l=1; l<L; ++l, ++X, ++Y) { *Y = *(Y-1) * *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                *Y++ = *X++;
                for (size_t l=1; l<L; ++l, ++X, ++Y) { *Y = *(Y-1) * *X; }
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *X; }
            for (size_t l=1; l<L; ++l)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *(Y-V) * *X; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, Y-=K*L-1)
                {
                    *Y = *X; X += K; Y += K;
                    for (size_t l=1; l<L; ++l, X+=K, Y+=K) { *Y = *(Y-K) * *X; }
                }
            }
        }
    }

    return 0;
}


int cumprod_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cumprod_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y++ = *X++;
        for (size_t l=1; l<L; ++l, ++X, ++Y) { *Y = *(Y-1) * *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                *Y++ = *X++;
                for (size_t l=1; l<L; ++l, ++X, ++Y) { *Y = *(Y-1) * *X; }
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *X; }
            for (size_t l=1; l<L; ++l)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *(Y-V) * *X; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, Y-=K*L-1)
                {
                    *Y = *X; X += K; Y += K;
                    for (size_t l=1; l<L; ++l, X+=K, Y+=K) { *Y = *(Y-K) * *X; }
                }
            }
        }
    }

    return 0;
}


int cumprod_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cumprod_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y++ = *X++; *Y++ = *X++;
        for (size_t l=1; l<L; ++l, X+=2, Y+=2)
        {
            *Y = *(Y-2)**X - *(Y-1)**(X+1);
            *(Y+1) = *(Y-1)**X + *(Y-2)**(X+1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                *Y++ = *X++; *Y++ = *X++;
                for (size_t l=1; l<L; ++l, ++X, ++Y) { *Y = *(Y-2) + *X; ++Y; *Y = *(Y-2) + *++X; }
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            for (size_t l=1; l<L; ++l)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *(Y-2*V) + *X; ++Y; *Y = *(Y-2*V) + *++X; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1), Y+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, Y-=2*K*L-2)
                {
                    *Y = *X; *(Y+1) = *(X+1); X += 2*K; Y += 2*K;
                    for (size_t l=1; l<L; ++l, X+=2*K, Y+=2*K) { *Y = *(Y-2*K) + *X; *(Y+1) = *(Y-2*K+1) + *(X+1); }
                }
            }
        }
    }

    return 0;
}


int cumprod_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cumprod_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y++ = *X++; *Y++ = *X++;
        for (size_t l=1; l<L; ++l, ++Y) { *Y = *(Y-2) + *X++; ++Y; *Y = *(Y-2) + *X++; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                *Y++ = *X++; *Y++ = *X++;
                for (size_t l=1; l<L; ++l, ++X, ++Y) { *Y = *(Y-2) + *X; ++Y; *Y = *(Y-2) + *++X; }
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            for (size_t l=1; l<L; ++l)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *(Y-2*V) + *X; ++Y; *Y = *(Y-2*V) + *++X; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1), Y+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, Y-=2*K*L-2)
                {
                    *Y = *X; *(Y+1) = *(X+1); X += 2*K; Y += 2*K;
                    for (size_t l=1; l<L; ++l, X+=2*K, Y+=2*K) { *Y = *(Y-2*K) + *X; *(Y+1) = *(Y-2*K+1) + *(X+1); }
                }
            }
        }
    }

    return 0;
}


int cumprod_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cumprod_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0 || L==1) {}
    else if (L==N)
    {
        ++X;
        for (size_t l=1; l<L; ++l, ++X) { *X *= *(X-1); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                ++X;
                for (size_t l=1; l<L; ++l, ++X) { *X *= *(X-1); }
            }
        }
        else if (G==1)
        {
            X += V;
            for (size_t l=1; l<L; ++l)
            {
                for (size_t v=0; v<V; ++v, ++X) { *X *= *(X-V); }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1)
                {
                    X += K;
                    for (size_t l=1; l<L; ++l, X+=K) { *X *= *(X-K); }
                }
            }
        }
    }

    return 0;
}


int cumprod_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cumprod_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0 || L==1) {}
    else if (L==N)
    {
        ++X;
        for (size_t l=1; l<L; ++l, ++X) { *X *= *(X-1); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                ++X;
                for (size_t l=1; l<L; ++l, ++X) { *X *= *(X-1); }
            }
        }
        else if (G==1)
        {
            X += V;
            for (size_t l=1; l<L; ++l)
            {
                for (size_t v=0; v<V; ++v, ++X) { *X *= *(X-V); }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1)
                {
                    X += K;
                    for (size_t l=1; l<L; ++l, X+=K) { *X *= *(X-K); }
                }
            }
        }
    }

    return 0;
}


int cumprod_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cumprod_inplace_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float xr, xi, c, d;

    if (N==0 || L==1) {}
    else if (L==N)
    {
        c = *X++; d = *X++;
        for (size_t l=1; l<L; ++l, ++X)
        {
            xr = *X; xi = *(X+1);
            *X = xr*c - xi*d;
            *++X = xr*d + xi*c;
            c = *(X-1); d = *X;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                c = *X++; d = *X++;
                for (size_t l=1; l<L; ++l, ++X)
                {
                    xr = *X; xi = *(X+1);
                    *X = xr*c - xi*d;
                    *++X = xr*d + xi*c;
                    c = *(X-1); d = *X;
                }
            }
        }
        else if (G==1)
        {
            X += 2*V;
            for (size_t l=1; l<L; ++l)
            {
                for (size_t v=0; v<V; ++v, ++X)
                {
                    c = *(X-2*V); d = *(X-2*V+1);
                    xr = *X; xi = *(X+1);
                    *X = xr*c - xi*d;
                    *++X = xr*d + xi*c;
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2)
                {
                    c = *X; d = *++X; X += 2*K-1;
                    for (size_t l=1; l<L; ++l, X+=2*K-1)
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


int cumprod_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cumprod_inplace_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double xr, xi, c, d;

    if (N==0 || L==1) {}
    else if (L==N)
    {
        c = *X++; d = *X++;
        for (size_t l=1; l<L; ++l, ++X)
        {
            xr = *X; xi = *(X+1);
            *X = xr*c - xi*d;
            *++X = xr*d + xi*c;
            c = *(X-1); d = *X;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                c = *X++; d = *X++;
                for (size_t l=1; l<L; ++l, ++X)
                {
                    xr = *X; xi = *(X+1);
                    *X = xr*c - xi*d;
                    *++X = xr*d + xi*c;
                    c = *(X-1); d = *X;
                }
            }
        }
        else if (G==1)
        {
            X += 2*V;
            for (size_t l=1; l<L; ++l)
            {
                for (size_t v=0; v<V; ++v, ++X)
                {
                    c = *(X-2*V); d = *(X-2*V+1);
                    xr = *X; xi = *(X+1);
                    *X = xr*c - xi*d;
                    *++X = xr*d + xi*c;
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2)
                {
                    c = *X; d = *++X; X += 2*K-1;
                    for (size_t l=1; l<L; ++l, X+=2*K-1)
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
