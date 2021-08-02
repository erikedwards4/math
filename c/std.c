//Vec2scalar (reduction) operation.
//Gets standard deviation of each row or col of X according to dim.
//For complex case, output is real.

//The use of cblas_snrm2 was not faster, and required in-place processing to optimize.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int std_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
int std_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
int std_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
int std_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);


int std_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased)
{
    if (dim>3u) { fprintf(stderr,"error in std_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = 1.0f/(float)L, den2 = (biased) ? den : 1.0f/(float)(L-1u);

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        float x, mn = 0.0f, sm2 = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=0u; l<L; ++l) { x = *--X - mn; sm2 += x*x; }
        *Y = sqrtf(sm2*den2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float x, mn, sm2;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mn = sm2 = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                X -= L;
                mn *= den;
                for (size_t l=0u; l<L; ++l, ++X) { x = *X - mn; sm2 += x*x; }
                *Y = sqrtf(sm2*den2);
            }
        }
        else if (G==1u)
        {
            float x, *mn;
            if (!(mn=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in std_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0u; l<L; ++l, mn-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++mn) { *mn += *X; }
            }
            X -= N;
            for (size_t v=0u; v<V; ++v, ++mn) { *mn *= den; }
            mn -= V;
            for (size_t v=0u; v<V; ++v, ++X, ++mn, ++Y) { x = *X - *mn; *Y = x*x; }
            mn -= V; Y -= V;
            for (size_t l=1u; l<L; ++l, mn-=V, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++mn, ++Y) { x = *X - *mn; *Y += x*x; }
            }
            for (size_t v=0u; v<V; ++v, ++Y) { *Y = sqrtf(*Y*den2); }
            free(mn);
        }
        else
        {
            float x, mn, sm2;
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, ++X, ++Y)
                {
                    mn = sm2 = 0.0f;
                    for (size_t l=0u; l<L-1u; ++l, X+=K) { mn += *X; }
                    mn += *X; mn *= den;
                    for (size_t l=0u; l<L-1u; ++l, X-=K) { x = *X - mn; sm2 += x*x; }
                    x = *X - mn; sm2 += x*x;
                    *Y = sqrtf(sm2*den2);
                }
            }
        }
    }

    return 0;
}


int std_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased)
{
    if (dim>3u) { fprintf(stderr,"error in std_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = 1.0/(double)L, den2 = (biased) ? den : 1.0/(double)(L-1u);

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        double x, mn = 0.0, sm2 = 0.0;
        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=0u; l<L; ++l) { x = *--X - mn; sm2 += x*x; }
        *Y = sqrt(sm2*den2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double x, mn, sm2;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mn = sm2 = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                X -= L;
                mn *= den;
                for (size_t l=0u; l<L; ++l, ++X) { x = *X - mn; sm2 += x*x; }
                *Y = sqrt(sm2*den2);
            }
        }
        else if (G==1u)
        {
            double x, *mn;
            if (!(mn=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in std_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0u; l<L; ++l, mn-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++mn) { *mn += *X; }
            }
            X -= N;
            for (size_t v=0u; v<V; ++v, ++mn) { *mn *= den; }
            mn -= V;
            for (size_t v=0u; v<V; ++v, ++X, ++mn, ++Y) { x = *X - *mn; *Y = x*x; }
            mn -= V; Y -= V;
            for (size_t l=1u; l<L; ++l, mn-=V, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++mn, ++Y) { x = *X - *mn; *Y += x*x; }
            }
            for (size_t v=0u; v<V; ++v, ++Y) { *Y = sqrt(*Y*den2); }
            free(mn);
        }
        else
        {
            double x, mn, sm2;
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, ++X, ++Y)
                {
                    mn = sm2 = 0.0;
                    for (size_t l=0u; l<L-1u; ++l, X+=K) { mn += *X; }
                    mn += *X; mn *= den;
                    for (size_t l=0u; l<L-1u; ++l, X-=K) { x = *X - mn; sm2 += x*x; }
                    x = *X - mn; sm2 += x*x;
                    *Y = sqrt(sm2*den2);
                }
            }
        }
    }

    return 0;
}


int std_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased)
{
    if (dim>3u) { fprintf(stderr,"error in std_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = 1.0f/(float)L, den2 = (biased) ? den : 1.0f/(float)(L-1u);
    float xr, xi;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        float mnr = 0.0f, mni = 0.0f, sm2 = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
        mnr *= den; mni *= den;
        X -= 2u*L;
        for (size_t l=0u; l<L; ++l, ++X) { xr = *X - mnr; xi = *++X - mni; sm2 += xr*xr + xi*xi; }
        *Y = sqrtf(sm2*den2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float mnr, mni, sm2;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mnr = mni = sm2 = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
                mnr *= den; mni *= den;
                X -= 2u*L;
                for (size_t l=0u; l<L; ++l, ++X) { xr = *X - mnr; xi = *++X - mni; sm2 += xr*xr + xi*xi; }
                *Y = sqrtf(sm2*den2);
            }
        }
        else
        {
            float mnr, mni, sm2;
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, ++Y)
                {
                    mnr = mni = sm2 = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u) { mnr += *X; mni += *++X; }
                    mnr *= den; mni *= den;
                    X -= 2u*K*L;
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u) { xr = *X - mnr; xi = *++X - mni; sm2 += xr*xr + xi*xi; }
                    *Y = sqrtf(sm2*den2);
                }
            }
        }
    }
    
    return 0;
}


int std_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased)
{
    if (dim>3u) { fprintf(stderr,"error in std_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = 1.0/(double)L, den2 = (biased) ? den : 1.0/(double)(L-1u);
    double xr, xi;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        double mnr = 0.0, mni = 0.0, sm2 = 0.0;
        for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
        mnr *= den; mni *= den;
        X -= 2u*L;
        for (size_t l=0u; l<L; ++l, ++X) { xr = *X - mnr; xi = *++X - mni; sm2 += xr*xr + xi*xi; }
        *Y = sqrt(sm2*den2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double mnr, mni, sm2;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mnr = mni = sm2 = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
                mnr *= den; mni *= den;
                X -= 2u*L;
                for (size_t l=0u; l<L; ++l, ++X) { xr = *X - mnr; xi = *++X - mni; sm2 += xr*xr + xi*xi; }
                *Y = sqrt(sm2*den2);
            }
        }
        else
        {
            double mnr, mni, sm2;
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, ++Y)
                {
                    mnr = mni = sm2 = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u) { mnr += *X; mni += *++X; }
                    mnr *= den; mni *= den;
                    X -= 2u*K*L;
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u) { xr = *X - mnr; xi = *++X - mni; sm2 += xr*xr + xi*xi; }
                    *Y = sqrt(sm2*den2);
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
