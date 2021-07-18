//Vec2scalar (reduction) operation.
//Gets harmonic mean for each vector in X along dim.
//This is the reciprocal of the mean of reciprocals for each vector.

//The same definition is used for complex numbers (as in Octave),
//but there are some publications on this topic, so check into later if needed.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int harmean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int harmean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int harmean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int harmean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int harmean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in harmean_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        float sm = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X) { sm += 1.0f / *X; }
        *Y = (float)L / sm;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float sm;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                sm = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { sm += 1.0f / *X; }
                *Y = (float)L / sm;
            }
        }
        else if (G==1u)
        {
            for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y = 1.0f / *X; }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += 1.0f / *X; }
            }
            for (size_t v=0u; v<V; ++v, ++Y) { *Y = (float)L / *Y; }
        }
        else
        {
            float sm;
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=K) { sm += 1.0f / *X; }
                    *Y = (float)L / sm;
                }
            }
        }
    }

    return 0;
}


int harmean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in harmean_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        double sm = 0.0;
        for (size_t l=0u; l<L; ++l, ++X) { sm += 1.0 / *X; }
        *Y = (double)L / sm;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double sm;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                sm = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { sm += 1.0 / *X; }
                *Y = (double)L / sm;
            }
        }
        else if (G==1u)
        {
            for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y = 1.0 / *X; }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += 1.0 / *X; }
            }
            for (size_t v=0u; v<V; ++v, ++Y) { *Y = (double)L / *Y; }
        }
        else
        {
            double sm;
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=K) { sm += 1.0 / *X; }
                    *Y = (double)L / sm;
                }
            }
        }
    }

    return 0;
}


int harmean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in harmean_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float xr, xi, sq, yr, yi;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        yr = yi = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X)
        {
            xr = *X; xi = *++X;
            sq = xr*xr + xi*xi;
            yr += xr/sq; yi += xi/sq;
        }
        sq = yr*yr + yi*yi;
        *Y = (float)L*yr/sq; *++Y = (float)L*yi/sq;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                yr = yi = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X)
                {
                    xr = *X; xi = *++X;
                    sq = xr*xr + xi*xi;
                    yr += xr/sq; yi += xi/sq;
                }
                sq = yr*yr + yi*yi;
                *Y = (float)L*yr/sq; *++Y = (float)L*yi/sq;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, ++Y)
                {
                    yr = yi = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u)
                    {
                        xr = *X; xi = *++X;
                        sq = xr*xr + xi*xi;
                        yr += xr/sq; yi += xi/sq;
                    }
                    sq = yr*yr + yi*yi;
                    *Y = (float)L*yr/sq; *++Y = (float)L*yi/sq;
                }
            }
        }
    }

    return 0;
}


int harmean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in harmean_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double xr, xi, sq, yr, yi;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        yr = yi = 0.0;
        for (size_t l=0u; l<L; ++l, ++X)
        {
            xr = *X; xi = *++X;
            sq = xr*xr + xi*xi;
            yr += xr/sq; yi += xi/sq;
        }
        sq = yr*yr + yi*yi;
        *Y = (double)L*yr/sq; *++Y = (double)L*yi/sq;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                yr = yi = 0.0;
                for (size_t l=0u; l<L; ++l, ++X)
                {
                    xr = *X; xi = *++X;
                    sq = xr*xr + xi*xi;
                    yr += xr/sq; yi += xi/sq;
                }
                sq = yr*yr + yi*yi;
                *Y = (double)L*yr/sq; *++Y = (double)L*yi/sq;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, ++Y)
                {
                    yr = yi = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u)
                    {
                        xr = *X; xi = *++X;
                        sq = xr*xr + xi*xi;
                        yr += xr/sq; yi += xi/sq;
                    }
                    sq = yr*yr + yi*yi;
                    *Y = (double)L*yr/sq; *++Y = (double)L*yi/sq;
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
