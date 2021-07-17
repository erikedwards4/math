//Vec2vec operation.
//Scales each vector in X to a range of 0 to 1 along dim, i.e., min-max scaling.
//If m1, then scales to a range of -1 to 1.
//This operates in-place.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int range1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char m1);
int range1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char m1);


int range1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char m1)
{
    if (dim>3u) { fprintf(stderr,"error in range1_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in range1_s: L (vec length) must be > 1\n"); return 1; }
    float mn, mx, rng;

    if (N==0u) {}
    else if (L==N)
    {
        mn = mx = *X++;
        for (size_t l=1u; l<L; ++l, ++X)
        {
            if (*X<mn) { mn = *X; }
            else if (*X>mx) { mx = *X; }
        }
        rng = mx - mn;
        if (m1)
        {
            for (size_t l=0u; l<L; ++l) { --X; *X = 2.0f*(*X-mn)/rng - 1.0f; }
        }
        else
        {
            for (size_t l=0u; l<L; ++l) { --X; *X = (*X-mn)/rng; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                mn = mx = *X++;
                for (size_t l=1u; l<L; ++l, ++X)
                {
                    if (*X<mn) { mn = *X; }
                    else if (*X>mx) { mx = *X; }
                }
                rng = mx - mn; X -= L;
                if (m1)
                {
                    for (size_t l=0u; l<L; ++l, ++X) { *X = 2.0f*(*X-mn)/rng - 1.0f; }
                }
                else
                {
                    for (size_t l=0u; l<L; ++l, ++X) { *X = (*X-mn)/rng; }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, ++X)
                {
                    mn = mx = *X; X += K;
                    for (size_t l=1u; l<L; ++l, X+=K)
                    {
                        if (*X<mn) { mn = *X; }
                        else if (*X>mx) { mx = *X; }
                    }
                    rng = mx - mn;
                    if (m1)
                    {
                        for (size_t l=0u; l<L; ++l) { X-=K; *X = 2.0f*(*X-mn)/rng - 1.0f; }
                    }
                    else
                    {
                        for (size_t l=0u; l<L; ++l) { X-=K; *X = (*X-mn)/rng; }
                    }
                }
            }
        }
    }

    return 0;
}


int range1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char m1)
{
    if (dim>3u) { fprintf(stderr,"error in range1_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in range1_d: L (vec length) must be > 1\n"); return 1; }
    double mn, mx, rng;

    if (N==0u) {}
    else if (L==N)
    {
        mn = mx = *X++;
        for (size_t l=1u; l<L; ++l, ++X)
        {
            if (*X<mn) { mn = *X; }
            else if (*X>mx) { mx = *X; }
        }
        rng = mx - mn;
        if (m1)
        {
            for (size_t l=0u; l<L; ++l) { --X; *X = 2.0*(*X-mn)/rng - 1.0; }
        }
        else
        {
            for (size_t l=0u; l<L; ++l) { --X; *X = (*X-mn)/rng; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                mn = mx = *X++;
                for (size_t l=1u; l<L; ++l, ++X)
                {
                    if (*X<mn) { mn = *X; }
                    else if (*X>mx) { mx = *X; }
                }
                rng = mx - mn; X -= L;
                if (m1)
                {
                    for (size_t l=0u; l<L; ++l, ++X) { *X = 2.0*(*X-mn)/rng - 1.0; }
                }
                else
                {
                    for (size_t l=0u; l<L; ++l, ++X) { *X = (*X-mn)/rng; }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, ++X)
                {
                    mn = mx = *X; X += K;
                    for (size_t l=1u; l<L; ++l, X+=K)
                    {
                        if (*X<mn) { mn = *X; }
                        else if (*X>mx) { mx = *X; }
                    }
                    rng = mx - mn;
                    if (m1)
                    {
                        for (size_t l=0u; l<L; ++l) { X-=K; *X = 2.0*(*X-mn)/rng - 1.0; }
                    }
                    else
                    {
                        for (size_t l=0u; l<L; ++l) { X-=K; *X = (*X-mn)/rng; }
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
