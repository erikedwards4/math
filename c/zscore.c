//Z-scores each vector in X along dim.
//This operates in-place.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int zscore_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int zscore_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int zscore_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int zscore_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);


int zscore_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in zscore_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (L<2) { fprintf(stderr,"error in zscore_s: L (vec length) must be > 1\n"); return 1; }
    const float den = 1.0f/L, den2 = (biased) ? den : 1.0f/(L-1);
    float mn = 0.0f, sd = 0.0f;

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=0u; l<L; ++l) { *--X -= mn; sd += *X**X; }
        sd = sqrtf(sd*den2);
        for (size_t l=0u; l<L; ++l, ++X) { *X /= sd; }
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
                mn = sd = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                mn *= den;
                for (size_t l=0u; l<L; ++l) { *--X -= mn; sd += *X**X; }
                sd = sqrtf(sd*den2);
                for (size_t l=0u; l<L; ++l, ++X) { *X /= sd; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1)
                {
                    mn = sd = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=K) { mn += *X; }
                    mn *= den;
                    for (size_t l=0u; l<L; ++l) { X-=K; *X -= mn; sd += *X**X; }
                    sd = sqrtf(sd*den2);
                    for (size_t l=0u; l<L; ++l, X+=K) { *X /= sd; }
                }
            }
        }
    }

    return 0;
}


int zscore_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in zscore_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (L<2) { fprintf(stderr,"error in zscore_d: L (vec length) must be > 1\n"); return 1; }
    const double den = 1.0/L, den2 = (biased) ? den : 1.0/(L-1);
    double mn = 0.0, sd = 0.0;

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=0u; l<L; ++l) { *--X -= mn; sd += *X**X; }
        sd = sqrt(sd*den2);
        for (size_t l=0u; l<L; ++l, ++X) { *X /= sd; }
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
                mn = sd = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                mn *= den;
                for (size_t l=0u; l<L; ++l) { *--X -= mn; sd += *X**X; }
                sd = sqrt(sd*den2);
                for (size_t l=0u; l<L; ++l, ++X) { *X /= sd; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1)
                {
                    mn = sd = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=K) { mn += *X; }
                    mn *= den;
                    for (size_t l=0u; l<L; ++l) { X-=K; *X -= mn; sd += *X**X; }
                    sd = sqrt(sd*den2);
                    for (size_t l=0u; l<L; ++l, X+=K) { *X /= sd; }
                }
            }
        }
    }

    return 0;
}


int zscore_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in zscore_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (L<2) { fprintf(stderr,"error in zscore_c: L (vec length) must be > 1\n"); return 1; }
    const float den = 1.0f/L, den2 = (biased) ? den : 1.0f/(L-1);
    float mnr = 0.0f, mni = 0.0f, sd = 0.0f;

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
        mnr *= den; mni *= den;
        for (size_t l=0u; l<L; ++l) { *--X -= mni; sd += *X**X; *--X -= mnr; sd += *X**X; }
        sd = sqrtf(sd*den2);
        for (size_t l=0u; l<2*L; ++l, ++X) { *X /= sd; }
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
                mnr = mni = sd = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
                mnr *= den; mni *= den;
                for (size_t l=0u; l<L; ++l) { *--X -= mni; sd += *X**X; *--X -= mnr; sd += *X**X; }
                sd = sqrtf(sd*den2);
                for (size_t l=0u; l<2*L; ++l, ++X) { *X /= sd; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2)
                {
                    mnr = mni = sd = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=2*K-1) { mnr += *X; mni += *++X; }
                    mnr *= den; mni *= den;
                    for (size_t l=0u; l<L; ++l) { X-=2*K-1; *X -= mni; sd += *X**X; *--X -= mnr; sd += *X**X; }
                    sd = sqrtf(sd*den2);
                    for (size_t l=0u; l<L; ++l, X+=2*K-1) { *X /= sd; *++X /= sd; }
                }
            }
        }
    }

    return 0;
}


int zscore_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in zscore_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (L<2) { fprintf(stderr,"error in zscore_c: L (vec length) must be > 1\n"); return 1; }
    const double den = 1.0/L, den2 = (biased) ? den : 1.0/(L-1);
    double mnr = 0.0, mni = 0.0, sd = 0.0;

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
        mnr *= den; mni *= den;
        for (size_t l=0u; l<L; ++l) { *--X -= mni; sd += *X**X; *--X -= mnr; sd += *X**X; }
        sd = sqrt(sd*den2);
        for (size_t l=0u; l<2*L; ++l, ++X) { *X /= sd; }
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
                mnr = mni = sd = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
                mnr *= den; mni *= den;
                for (size_t l=0u; l<L; ++l) { *--X -= mni; sd += *X**X; *--X -= mnr; sd += *X**X; }
                sd = sqrt(sd*den2);
                for (size_t l=0u; l<2*L; ++l, ++X) { *X /= sd; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2)
                {
                    mnr = mni = sd = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=2*K-1) { mnr += *X; mni += *++X; }
                    mnr *= den; mni *= den;
                    for (size_t l=0u; l<L; ++l) { X-=2*K-1; *X -= mni; sd += *X**X; *--X -= mnr; sd += *X**X; }
                    sd = sqrt(sd*den2);
                    for (size_t l=0u; l<L; ++l, X+=2*K-1) { *X /= sd; *++X /= sd; }
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
