//Vec2scalar (reduction) operation.
//Gets mean for each vector in X along dim.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int mean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in mean_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = 1.0f / (float)L;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (L<100u)
        {
            *Y = 0.0f;
            for (size_t l=L; l>0u; --l, ++X) { *Y += *X; }
            *Y *= den;
        }
        else
        {
            float sm = 0.0f;
            for (size_t l=L; l>0u; --l, ++X) { sm += *X; }
            *Y = sm * den;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (L<10u)
            {
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    *Y = 0.0f;
                    for (size_t l=L; l>0u; --l, ++X) { *Y += *X; }
                    *Y *= den;
                }
            }
            else
            {
                float sm;
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=L; l>0u; --l, ++X) { sm += *X; }
                    *Y = sm * den;
                }
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += *X; }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y *= den; }
        }
        else
        {
            float sm;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=K) { sm += *X; }
                    *Y = sm * den;
                }
            }
        }
    }

    return 0;
}


int mean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in mean_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = 1.0 / (double)L;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (L<100u)
        {
            *Y = 0.0;
            for (size_t l=L; l>0u; --l, ++X) { *Y += *X; }
            *Y *= den;
        }
        else
        {
            double sm = 0.0;
            for (size_t l=L; l>0u; --l, ++X) { sm += *X; }
            *Y = sm * den;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (L<10u)
            {
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    *Y = 0.0;
                    for (size_t l=L; l>0u; --l, ++X) { *Y += *X; }
                    *Y *= den;
                }
            }
            else
            {
                double sm;
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=L; l>0u; --l, ++X) { sm += *X; }
                    *Y = sm * den;
                }
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += *X; }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y *= den; }
        }
        else
        {
            double sm;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=L; l>0u; --l, X+=K) { sm += *X; }
                    *Y = sm * den;
                }
            }
        }
    }

    return 0;
}


int mean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in mean_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = 1.0f / (float)L;
    float yr, yi;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        yr = yi = 0.0f;
        for (size_t l=L; l>0u; --l, ++X) { yr += *X; yi += *++X; }
        *Y = yr * den; *++Y = yi * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++Y)
            {
                yr = yi = 0.0f;
                for (size_t l=L; l>0u; --l, ++X) { yr += *X; yi += *++X; }
                *Y = yr * den; *++Y = yi * den;
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            Y -= 2u*V;
            for (size_t l=1u; l<L; ++l, Y-=2u*V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += *X; *++Y += *++X; }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y *= den; *++Y *= den; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    yr = yi = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u) { yr += *X; yi += *++X; }
                    *Y = yr * den; *++Y = yi * den;
                }
            }
        }
    }

    return 0;
}


int mean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in mean_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = 1.0 / (double)L;
    double yr, yi;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        yr = yi = 0.0;
        for (size_t l=L; l>0u; --l, ++X) { yr += *X; yi += *++X; }
        *Y = yr * den; *++Y = yi * den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++Y)
            {
                yr = yi = 0.0;
                for (size_t l=L; l>0u; --l, ++X) { yr += *X; yi += *++X; }
                *Y = yr * den; *++Y = yi * den;
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            Y -= 2u*V;
            for (size_t l=1u; l<L; ++l, Y-=2u*V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += *X; *++Y += *++X; }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y *= den; *++Y *= den; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    yr = yi = 0.0;
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u) { yr += *X; yi += *++X; }
                    *Y = yr * den; *++Y = yi * den;
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
