//"G"-scores each vector in X along dim.
//For each vector, z-scores the log values and then takes exp.
//This operates in-place.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int gscore_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int gscore_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int gscore_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in gscore_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (L<2) { fprintf(stderr,"error in gscore_s: L (vec length) must be > 1\n"); return 1; }
    const float den = 1.0f/L, den2 = den;
    float mn = 0.0f, sd = 0.0f;

    if (N==0) {}
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X) { *X = logf(*X); mn += *X; }
        mn *= den;
        for (size_t l=0; l<L; ++l) { *--X -= mn; sd += *X**X; }
        sd = sqrtf(sd*den2);
        for (size_t l=0; l<L; ++l, ++X) { *X = expf(*X/sd); }
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
                mn = sd = 0.0f;
                for (size_t l=0; l<L; ++l, ++X) { *X = logf(*X); mn += *X; }
                mn *= den;
                for (size_t l=0; l<L; ++l) { *--X -= mn; sd += *X**X; }
                sd = sqrtf(sd*den2);
                for (size_t l=0; l<L; ++l, ++X) { *X = expf(*X/sd); }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1)
                {
                    mn = sd = 0.0f;
                    for (size_t l=0; l<L; ++l, X+=K) { *X = logf(*X); mn += *X; }
                    mn *= den;
                    for (size_t l=0; l<L; ++l) { X-=K; *X -= mn; sd += *X**X; }
                    sd = sqrtf(sd*den2);
                    for (size_t l=0; l<L; ++l, X+=K) { *X = expf(*X/sd); }
                }
            }
        }
    }

    return 0;
}


int gscore_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in gscore_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (L<2) { fprintf(stderr,"error in gscore_s: L (vec length) must be > 1\n"); return 1; }
    const double den = 1.0/L, den2 = den;
    double mn = 0.0, sd = 0.0;

    if (N==0) {}
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X) { *X = log(*X); mn += *X; }
        mn *= den;
        for (size_t l=0; l<L; ++l) { *--X -= mn; sd += *X**X; }
        sd = sqrt(sd*den2);
        for (size_t l=0; l<L; ++l, ++X) { *X = exp(*X/sd); }
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
                mn = sd = 0.0;
                for (size_t l=0; l<L; ++l, ++X) { *X = log(*X); mn += *X; }
                mn *= den;
                for (size_t l=0; l<L; ++l) { *--X -= mn; sd += *X**X; }
                sd = sqrt(sd*den2);
                for (size_t l=0; l<L; ++l, ++X) { *X = exp(*X/sd); }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1)
                {
                    mn = sd = 0.0;
                    for (size_t l=0; l<L; ++l, X+=K) { *X = log(*X); mn += *X; }
                    mn *= den;
                    for (size_t l=0; l<L; ++l) { X-=K; *X -= mn; sd += *X**X; }
                    sd = sqrt(sd*den2);
                    for (size_t l=0; l<L; ++l, X+=K) { *X = exp(*X/sd); }
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
