//Zeros the mean of each vector in X along dim.
//For each vector, estimates the mean and then subtracts it from each element.
//This operates in-place.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mean0_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mean0_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mean0_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mean0_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int mean0_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean0_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / L;
    float mn = 0.0f;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=0u; l<L; ++l) { *--X -= mn; }
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
                for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                mn *= den; X -= L;
                for (size_t l=0u; l<L; ++l, ++X) { *X -= mn; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, ++X)
                {
                    mn = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=K) { mn += *X; }
                    mn *= den;
                    for (size_t l=0u; l<L; ++l) { X-=K; *X -= mn; }
                }
            }
        }
    }

    return 0;
}


int mean0_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean0_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / L;
    double mn = 0.0;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=0u; l<L; ++l) { *--X -= mn; }
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
                for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                mn *= den; X -= L;
                for (size_t l=0u; l<L; ++l, ++X) { *X -= mn; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, ++X)
                {
                    mn = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=K) { mn += *X; }
                    mn *= den;
                    for (size_t l=0u; l<L; ++l) { X-=K; *X -= mn; }
                }
            }
        }
    }

    return 0;
}


int mean0_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean0_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / L;
    float mnr = 0.0f, mni = 0.0f;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2*N; ++n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
        mnr *= den; mni *= den;
        for (size_t l=0u; l<L; ++l) { *--X -= mnr; *--X -= mni; }
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
                mnr = mni = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
                mnr *= den; mni *= den; X -= 2*L;
                for (size_t l=0u; l<L; ++l, ++X) { *X -= mnr; *++X -= mni; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X+=2)
                {
                    mnr = mni = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=2*K-1) { mnr += *X; mni += *++X; }
                    mnr *= den; mni *= den;
                    for (size_t l=0u; l<L; ++l) { X-=2*K-1; *X -= mni; *--X -= mnr; }
                }
            }
        }
    }

    return 0;
}


int mean0_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean0_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / L;
    double mnr = 0.0, mni = 0.0;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2*N; ++n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
        mnr *= den; mni *= den;
        for (size_t l=0u; l<L; ++l) { *--X -= mnr; *--X -= mni; }
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
                mnr = mni = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
                mnr *= den; mni *= den; X -= 2*L;
                for (size_t l=0u; l<L; ++l, ++X) { *X -= mnr; *++X -= mni; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X+=2)
                {
                    mnr = mni = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=2*K-1) { mnr += *X; mni += *++X; }
                    mnr *= den; mni *= den;
                    for (size_t l=0u; l<L; ++l) { X-=2*K-1; *X -= mni; *--X -= mnr; }
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
