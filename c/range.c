//Vec2scalar (reduction) operation.
//Gets range of each vector in X along dim.
//The range is the 100th minus the 0th percentile (max - min).

//There is no need for an in-place version here.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int range_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int range_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int range_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in range_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float mn, mx;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        mn = mx = *X++;
        for (size_t l=1u; l<L; ++l, ++X)
        {
            if (*X<mn) { mn = *X; }
            else if (*X>mx) { mx = *X; }
        }
        *Y = mx - mn;
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
                mn = mx = *X++;
                for (size_t l=1u; l<L; ++l, ++X)
                {
                    if (*X<mn) { mn = *X; }
                    else if (*X>mx) { mx = *X; }
                }
                *Y = mx - mn;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    mn = mx = *X; X += K;
                    for (size_t l=1u; l<L; ++l, X+=K)
                    {
                        if (*X<mn) { mn = *X; }
                        else if (*X>mx) { mx = *X; }
                    }
                    *Y = mx - mn;
                }
            }
        }
    }

    return 0;
}


int range_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in range_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double mn, mx;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        mn = mx = *X++;
        for (size_t l=1u; l<L; ++l, ++X)
        {
            if (*X<mn) { mn = *X; }
            else if (*X>mx) { mx = *X; }
        }
        *Y = mx - mn;
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
                mn = mx = *X++;
                for (size_t l=1u; l<L; ++l, ++X)
                {
                    if (*X<mn) { mn = *X; }
                    else if (*X>mx) { mx = *X; }
                }
                *Y = mx - mn;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    mn = mx = *X; X += K;
                    for (size_t l=1u; l<L; ++l, X+=K)
                    {
                        if (*X<mn) { mn = *X; }
                        else if (*X>mx) { mx = *X; }
                    }
                    *Y = mx - mn;
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
