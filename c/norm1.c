//Vec2scalar (reduction) operation.
//Gets L1 norm (taxicab) for each vector in X along dim.
//This is the sum of absolute-values for each vector in X.

//For complex case, output is real. This does not use cblas_scasum or cblas_dzasum,
//because those sum |xr|+|xi|, which is not the definition of the L1 norm.
//Thus, this code is identical to asum for real-valued case, but not for complex.

//I could not confirm any case where cblas_?asum is faster than direct sum.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int norm1_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int norm1_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int norm1_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int norm1_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int norm1_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in norm1_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = fabsf(*X); }
    }
    else if (L==N)
    {
        float sm = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X) { sm += fabsf(*X); }
        *Y = sm;
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
                for (size_t l=0u; l<L; ++l, ++X) { sm += fabsf(*X); }
                *Y = sm;
            }
        }
        else if (G==1u)
        {
            for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += fabsf(*X); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += fabsf(*X); }
            }
        }
        else
        {
            float sm;
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=K) { sm += fabsf(*X); }
                    *Y = sm;
                }
            }
        }
    }

    return 0;
}


int norm1_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in norm1_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = fabs(*X); }
    }
    else if (L==N)
    {
        double sm = 0.0;
        for (size_t l=0u; l<L; ++l, ++X) { sm += fabs(*X); }
        *Y = sm;
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
                for (size_t l=0u; l<L; ++l, ++X) { sm += fabs(*X); }
                *Y = sm;
            }
        }
        else if (G==1u)
        {
            for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += fabs(*X); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += fabs(*X); }
            }
        }
        else
        {
            double sm;
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=K) { sm += fabs(*X); }
                    *Y = sm;
                }
            }
        }
    }

    return 0;
}


int norm1_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in norm1_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, X+=2, ++Y) { *Y = sqrtf(*X**X + *(X+1)**(X+1)); }
    }
    else if (L==N)
    {
        float sm = 0.0f;
        for (size_t l=0u; l<L; ++l, X+=2) { sm += sqrtf(*X**X + *(X+1)**(X+1)); }
        *Y = sm;
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
                for (size_t l=0u; l<L; ++l, X+=2) { sm += sqrtf(*X**X + *(X+1)**(X+1)); }
                *Y = sm;
            }
        }
        else if (G==1u)
        {
            for (size_t v=0u; v<V; ++v, X+=2, ++Y) { *Y = sqrtf(*X**X + *(X+1)**(X+1)); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, X+=2, ++Y) { *Y += sqrtf(*X**X + *(X+1)**(X+1)); }
            }
        }
        else
        {
            float sm;
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=2u*K) { sm += sqrtf(*X**X + *(X+1)**(X+1)); }
                    *Y = sm;
                }
            }
        }
    }

    return 0;
}


int norm1_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in norm1_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, X+=2, ++Y) { *Y = sqrt(*X**X + *(X+1)**(X+1)); }
    }
    else if (L==N)
    {
        double sm = 0.0;
        for (size_t l=0u; l<L; ++l, X+=2) { sm += sqrt(*X**X + *(X+1)**(X+1)); }
        *Y = sm;
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
                for (size_t l=0u; l<L; ++l, X+=2) { sm += sqrt(*X**X + *(X+1)**(X+1)); }
                *Y = sm;
            }
        }
        else if (G==1u)
        {
            for (size_t v=0u; v<V; ++v, X+=2, ++Y) { *Y = sqrt(*X**X + *(X+1)**(X+1)); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, X+=2, ++Y) { *Y += sqrt(*X**X + *(X+1)**(X+1)); }
            }
        }
        else
        {
            double sm;
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=2u*K) { sm += sqrt(*X**X + *(X+1)**(X+1)); }
                    *Y = sm;
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
