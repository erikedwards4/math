//Centers each vector in X along dim by subtracting the mean of log values,
//and then taking the exp, resulting in a geometric mean of 1.
//This is appropriate for all-positive X.
//This operates in-place.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int geomean1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int geomean1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int geomean1_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int geomean1_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int geomean1_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in geomean1_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / L;
    float mn = 0.0f;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *X = 1.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { *X = logf(*X); mn += *X; }
        mn *= den;
        for (size_t l=0u; l<L; ++l) { --X; *X = expf(*X-mn); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                mn = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { *X = logf(*X); mn += *X; }
                mn *= den; X -= L;
                for (size_t l=0u; l<L; ++l, ++X) { *X = expf(*X-mn); }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, ++X)
                {
                    mn = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=K) { *X = logf(*X); mn += *X; }
                    mn *= den;
                    for (size_t l=0u; l<L; ++l) { X-=K; *X = expf(*X-mn); }
                }
            }
        }
    }

    return 0;
}


int geomean1_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in geomean1_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / L;
    double mn = 0.0;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *X = 1.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { *X = log(*X); mn += *X; }
        mn *= den;
        for (size_t l=0u; l<L; ++l) { --X; *X = exp(*X-mn); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                mn = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { *X = log(*X); mn += *X; }
                mn *= den; X -= L;
                for (size_t l=0u; l<L; ++l, ++X) { *X = exp(*X-mn); }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, ++X)
                {
                    mn = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=K) { *X = log(*X); mn += *X; }
                    mn *= den;
                    for (size_t l=0u; l<L; ++l) { X-=K; *X = exp(*X-mn); }
                }
            }
        }
    }

    return 0;
}


int geomean1_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in geomean1_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / L;
    _Complex float x, mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *X = 1.0f; *++X = 0.0f; }
    }
    else if (L==N)
    {
        mn = 0.0f + 0.0if;
        for (size_t l=0u; l<L; ++l, ++X)
        {
            x = clogf(*X + 1.0if**(X+1));
            mn += x;
            *X = *(float *)&x; *++X = *((float *)&x+1);
        }
        mn *= den; X -= 2*L;
        for (size_t l=0u; l<L; ++l, ++X)
        {
            x = *X + 1.0if**(X+1);
            x = cexpf(x-mn);
            *X = *(float *)&x; *++X = *((float *)&x+1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                mn = 0.0f + 0.0if;
                for (size_t l=0u; l<L; ++l, ++X)
                {
                    x = clogf(*X + 1.0if**(X+1));
                    mn += x;
                    *X = *(float *)&x; *++X = *((float *)&x+1);
                }
                mn *= den; X -= 2*L;
                for (size_t l=0u; l<L; ++l, ++X)
                {
                    x = *X + 1.0if**(X+1);
                    x = cexpf(x-mn);
                    *X = *(float *)&x; *++X = *((float *)&x+1);
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2)
                {
                    mn = 0.0f + 0.0if;
                    for (size_t l=0u; l<L; ++l, X+=2*K-1)
                    {
                        x = clogf(*X + 1.0if**(X+1));
                        mn += x;
                        *X = *(float *)&x; *++X = *((float *)&x+1);
                    }
                    mn *= den; X -= 2*K*L;
                    for (size_t l=0u; l<L; ++l, X+=2*K-1)
                    {
                        x = *X + 1.0if**(X+1);
                        x = cexpf(x-mn);
                        *X = *(float *)&x; *++X = *((float *)&x+1);
                    }
                }
            }
        }
    }

    return 0;
}


int geomean1_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in geomean1_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / L;
    _Complex double x, mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *X = 1.0; *++X = 0.0; }
    }
    else if (L==N)
    {
        mn = 0.0 + 0.0i;
        for (size_t l=0u; l<L; ++l, ++X)
        {
            x = clog(*X + 1.0i**(X+1));
            mn += x;
            *X = *(double *)&x; *++X = *((double *)&x+1);
        }
        mn *= den; X -= 2*L;
        for (size_t l=0u; l<L; ++l, ++X)
        {
            x = *X + 1.0i**(X+1);
            x = cexp(x-mn);
            *X = *(double *)&x; *++X = *((double *)&x+1);
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                mn = 0.0 + 0.0i;
                for (size_t l=0u; l<L; ++l, ++X)
                {
                    x = clog(*X + 1.0i**(X+1));
                    mn += x;
                    *X = *(double *)&x; *++X = *((double *)&x+1);
                }
                mn *= den; X -= 2*L;
                for (size_t l=0u; l<L; ++l, ++X)
                {
                    x = *X + 1.0i**(X+1);
                    x = cexp(x-mn);
                    *X = *(double *)&x; *++X = *((double *)&x+1);
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2)
                {
                    mn = 0.0 + 0.0i;
                    for (size_t l=0u; l<L; ++l, X+=2*K-1)
                    {
                        x = clog(*X + 1.0i**(X+1));
                        mn += x;
                        *X = *(double *)&x; *++X = *((double *)&x+1);
                    }
                    mn *= den; X -= 2*K*L;
                    for (size_t l=0u; l<L; ++l, X+=2*K-1)
                    {
                        x = *X + 1.0i**(X+1);
                        x = cexp(x-mn);
                        *X = *(double *)&x; *++X = *((double *)&x+1);
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
