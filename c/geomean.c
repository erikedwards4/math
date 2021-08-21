//Vec2scalar (reduction) operation.
//Gets geometric mean for each vector in X along dim.
//This is the Lth root of the prod for each vector.
//This is also the exp of the mean of logs for each vector.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int geomean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int geomean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int geomean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int geomean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int geomean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in geomean_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = 1.0f / (float)L;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        float sm = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X) { sm += logf(*X); }
        *Y = expf(sm*den);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float sm;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                sm = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { sm += logf(*X); }
                *Y = expf(sm*den);
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = logf(*X); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += logf(*X); }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y = expf(*Y*den); }
        }
        else
        {
            float sm;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=K) { sm += logf(*X); }
                    *Y = expf(sm*den);
                }
            }
        }
    }

    return 0;
}


int geomean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in geomean_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = 1.0 / (double)L;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        double sm = 0.0;
        for (size_t l=0u; l<L; ++l, ++X) { sm += log(*X); }
        *Y = exp(sm*den);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double sm;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                sm = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { sm += log(*X); }
                *Y = exp(sm*den);
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = log(*X); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += log(*X); }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y = exp(*Y*den); }
        }
        else
        {
            double sm;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=K) { sm += log(*X); }
                    *Y = exp(sm*den);
                }
            }
        }
    }

    return 0;
}


int geomean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in geomean_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = 1.0f / (float)L;
    _Complex float y;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        y = 0.0f + 0.0if;
        for (size_t l=0u; l<L; ++l, X+=2)
        {
            y += clogf(*X + 1.0if**(X+1));
        }
        y = cexpf(y*den);
        *Y = *(float *)&y; *++Y = *((float *)&y+1);
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
                y = 0.0f + 0.0if;
                for (size_t l=0u; l<L; ++l, X+=2)
                {
                    y += clogf(*X + 1.0if**(X+1));
                }
                y = cexpf(y*den);
                *Y = *(float *)&y; *++Y = *((float *)&y+1);
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u)
                {
                    y = 0.0f + 0.0if;
                    for (size_t l=0u; l<L; ++l, X+=2u*K, ++Y)
                    {
                        y += clogf(*X + 1.0if**(X+1));
                    }
                    y = cexpf(y*den);
                    *Y = *(float *)&y; *++Y = *((float *)&y+1);
                }
            }
        }
    }

    return 0;
}


int geomean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in geomean_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = 1.0 / (double)L;
    _Complex double y;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        y = 0.0 + 0.0i;
        for (size_t l=0u; l<L; ++l, X+=2)
        {
            y += clog(*X + 1.0i**(X+1));
        }
        y = cexp(y*den);
        *Y = *(double *)&y; *++Y = *((double *)&y+1);
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
                y = 0.0 + 0.0i;
                for (size_t l=0u; l<L; ++l, X+=2)
                {
                    y += clog(*X + 1.0i**(X+1));
                }
                y = cexp(y*den);
                *Y = *(double *)&y; *++Y = *((double *)&y+1);
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u)
                {
                    y = 0.0 + 0.0i;
                    for (size_t l=0u; l<L; ++l, X+=2u*K, ++Y)
                    {
                        y += clog(*X + 1.0i**(X+1));
                    }
                    y = cexp(y*den);
                    *Y = *(double *)&y; *++Y = *((double *)&y+1);
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
