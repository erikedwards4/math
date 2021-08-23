//Zeros the mean of each vector in X along dim.
//For each vector, estimates the mean and then subtracts it from each element.
//This operates in-place.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mean0_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int mean0_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int mean0_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int mean0_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int mean0_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in mean0_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float mn = 0.0f;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X) { mn += *X; }
        mn /= (float)L;
        for (size_t l=L; l>0u; --l) { *--X -= mn; }
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
                mn = 0.0f;
                for (size_t l=L; l>0u; --l, ++X) { mn += *X; }
                mn /= (float)L; X -= L;
                for (size_t l=L; l>0u; --l, ++X) { *X -= mn; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    mn = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=K) { mn += *X; }
                    mn /= (float)L;
                    for (size_t l=L; l>0u; --l) { X-=K; *X -= mn; }
                }
            }
        }
    }

    return 0;
}


int mean0_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in mean0_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double mn = 0.0;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X) { mn += *X; }
        mn /= (double)L;
        for (size_t l=L; l>0u; --l) { *--X -= mn; }
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
                mn = 0.0;
                for (size_t l=L; l>0u; --l, ++X) { mn += *X; }
                mn /= (double)L; X -= L;
                for (size_t l=L; l>0u; --l, ++X) { *X -= mn; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    mn = 0.0;
                    for (size_t l=L; l>0u; --l, X+=K) { mn += *X; }
                    mn /= (double)L;
                    for (size_t l=L; l>0u; --l) { X-=K; *X -= mn; }
                }
            }
        }
    }

    return 0;
}


int mean0_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in mean0_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float mnr = 0.0f, mni = 0.0f;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X) { mnr += *X; mni += *++X; }
        mnr /= (float)L; mni /= (float)L;
        for (size_t l=L; l>0u; --l) { *--X -= mnr; *--X -= mni; }
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
                mnr = mni = 0.0f;
                for (size_t l=L; l>0u; --l, ++X) { mnr += *X; mni += *++X; }
                mnr /= (float)L; mni /= (float)L; X -= 2u*L;
                for (size_t l=L; l>0u; --l, ++X) { *X -= mnr; *++X -= mni; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X+=2)
                {
                    mnr = mni = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u) { mnr += *X; mni += *++X; }
                    mnr /= (float)L; mni /= (float)L;
                    for (size_t l=L; l>0u; --l) { X-=2u*K-1u; *X -= mni; *--X -= mnr; }
                }
            }
        }
    }

    return 0;
}


int mean0_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in mean0_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double mnr = 0.0, mni = 0.0;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X) { mnr += *X; mni += *++X; }
        mnr /= (double)L; mni /= (double)L;
        for (size_t l=L; l>0u; --l) { *--X -= mnr; *--X -= mni; }
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
                mnr = mni = 0.0;
                for (size_t l=L; l>0u; --l, ++X) { mnr += *X; mni += *++X; }
                mnr /= (double)L; mni /= (double)L; X -= 2u*L;
                for (size_t l=L; l>0u; --l, ++X) { *X -= mnr; *++X -= mni; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X+=2)
                {
                    mnr = mni = 0.0;
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u) { mnr += *X; mni += *++X; }
                    mnr /= (double)L; mni /= (double)L;
                    for (size_t l=L; l>0u; --l) { X-=2u*K-1u; *X -= mni; *--X -= mnr; }
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
