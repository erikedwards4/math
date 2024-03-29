//"G"-scores each vector in X along dim.
//For each vector, z-scores the log values and then takes exp.

#include <stdio.h>
#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int gscore_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in gscore_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in gscore_s: L (vec length) must be > 1\n"); return 1; }
    const float den = 1.0f/(float)L, den2 = den;
    float mn = 0.0f, sd = 0.0f;

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = logf(*X); mn += *Y; }
        mn *= den;
        for (size_t l=L; l>0u; --l) { *--Y -= mn; sd += *Y**Y; }
        sd = sqrtf(sd*den2);
        for (size_t l=L; l>0u; --l, ++Y) { *Y = expf(*Y/sd); }
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
                mn = sd = 0.0f;
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = logf(*X); mn += *Y; }
                mn *= den;
                for (size_t l=L; l>0u; --l) { *--Y -= mn; sd += *Y**Y; }
                sd = sqrtf(sd*den2);
                for (size_t l=L; l>0u; --l, ++Y) { *Y = expf(*Y/sd); }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                {
                    mn = sd = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=K, Y+=K) { *Y = logf(*X); mn += *Y; }
                    mn *= den;
                    for (size_t l=L; l>0u; --l) { Y-=K; *Y -= mn; sd += *Y**Y; }
                    sd = sqrtf(sd*den2);
                    for (size_t l=L; l>0u; --l, Y+=K) { *Y = expf(*Y/sd); }
                }
            }
        }
    }

    return 0;
}


int gscore_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in gscore_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in gscore_d: L (vec length) must be > 1\n"); return 1; }
    const double den = 1.0/(double)L, den2 = den;
    double mn = 0.0, sd = 0.0;

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = log(*X); mn += *Y; }
        mn *= den;
        for (size_t l=L; l>0u; --l) { *--Y -= mn; sd += *Y**Y; }
        sd = sqrt(sd*den2);
        for (size_t l=L; l>0u; --l, ++Y) { *Y = exp(*Y/sd); }
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
                mn = sd = 0.0;
                for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = log(*X); mn += *Y; }
                mn *= den;
                for (size_t l=L; l>0u; --l) { *--Y -= mn; sd += *Y**Y; }
                sd = sqrt(sd*den2);
                for (size_t l=L; l>0u; --l, ++Y) { *Y = exp(*Y/sd); }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                {
                    mn = sd = 0.0;
                    for (size_t l=L; l>0u; --l, X+=K, Y+=K) { *Y = log(*X); mn += *Y; }
                    mn *= den;
                    for (size_t l=L; l>0u; --l) { Y-=K; *Y -= mn; sd += *Y**Y; }
                    sd = sqrt(sd*den2);
                    for (size_t l=L; l>0u; --l, Y+=K) { *Y = exp(*Y/sd); }
                }
            }
        }
    }

    return 0;
}


int gscore_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in gscore_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in gscore_inplace_s: L (vec length) must be > 1\n"); return 1; }
    const float den = 1.0f/(float)L, den2 = den;
    float mn = 0.0f, sd = 0.0f;

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X) { *X = logf(*X); mn += *X; }
        mn *= den;
        for (size_t l=L; l>0u; --l) { *--X -= mn; sd += *X**X; }
        sd = sqrtf(sd*den2);
        for (size_t l=L; l>0u; --l, ++X) { *X = expf(*X/sd); }
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
                mn = sd = 0.0f;
                for (size_t l=L; l>0u; --l, ++X) { *X = logf(*X); mn += *X; }
                mn *= den;
                for (size_t l=L; l>0u; --l) { *--X -= mn; sd += *X**X; }
                sd = sqrtf(sd*den2);
                for (size_t l=L; l>0u; --l, ++X) { *X = expf(*X/sd); }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u)
                {
                    mn = sd = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=K) { *X = logf(*X); mn += *X; }
                    mn *= den;
                    for (size_t l=L; l>0u; --l) { X-=K; *X -= mn; sd += *X**X; }
                    sd = sqrtf(sd*den2);
                    for (size_t l=L; l>0u; --l, X+=K) { *X = expf(*X/sd); }
                }
            }
        }
    }

    return 0;
}


int gscore_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in gscore_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in gscore_inplace_s: L (vec length) must be > 1\n"); return 1; }
    const double den = 1.0/(double)L, den2 = den;
    double mn = 0.0, sd = 0.0;

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X) { *X = log(*X); mn += *X; }
        mn *= den;
        for (size_t l=L; l>0u; --l) { *--X -= mn; sd += *X**X; }
        sd = sqrt(sd*den2);
        for (size_t l=L; l>0u; --l, ++X) { *X = exp(*X/sd); }
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
                mn = sd = 0.0;
                for (size_t l=L; l>0u; --l, ++X) { *X = log(*X); mn += *X; }
                mn *= den;
                for (size_t l=L; l>0u; --l) { *--X -= mn; sd += *X**X; }
                sd = sqrt(sd*den2);
                for (size_t l=L; l>0u; --l, ++X) { *X = exp(*X/sd); }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u)
                {
                    mn = sd = 0.0;
                    for (size_t l=L; l>0u; --l, X+=K) { *X = log(*X); mn += *X; }
                    mn *= den;
                    for (size_t l=L; l>0u; --l) { X-=K; *X -= mn; sd += *X**X; }
                    sd = sqrt(sd*den2);
                    for (size_t l=L; l>0u; --l, X+=K) { *X = exp(*X/sd); }
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
