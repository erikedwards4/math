//Vec2scalar (reduction) operation.
//Gets skewness for each vector in X along dim.

//For complex case, output is complex.
//I follow the Octave convention for complex skewness (but see literature for other ideas later).

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int skewness_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int skewness_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int skewness_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int skewness_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);


int skewness_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in skewness_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / L;
    const float w = (biased) ? sqrtf(L) : L*sqrtf(L-1)/(L-2);

    if (N==0u) {}
    else if (L<3) { fprintf(stderr,"error in skewness_s: L must be > 2\n"); return 1; }
    else if (L==N)
    {
        float x, x2, mn = 0.0f, sm2 = 0.0f, sm3 = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=0u; l<L; ++l) { x = *--X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; }
        *Y = w * sm3 / (sm2*sqrtf(sm2));
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            float x, x2, mn, sm2, sm3;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mn = sm2 = sm3 = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                X -= L;
                mn *= den;
                for (size_t l=0u; l<L; ++l, ++X) { x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; }
                *Y = w * sm3 / (sm2*sqrtf(sm2));
            }
        }
        else if (G==1)
        {
            float x, x2, *sm2, *sm3;
            if (!(sm2=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm3=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += *X; }
            }
            X -= N;
            for (size_t v=0u; v<V; ++v, ++Y) { *Y *= den; }
            Y -= V;
            for (size_t l=0u; l<L; ++l, sm2-=V, sm3-=V, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++sm2, ++sm3, ++Y) { x = *X - *Y; x2 = x*x; *sm2 += x2; *sm3 += x*x2; }
            }
            for (size_t v=0u; v<V; ++v, ++sm2, ++sm3, ++Y) { *Y = w * *sm3 / (*sm2*sqrtf(*sm2)); }
            sm2 -= V; sm3 -= V;
            free(sm2); free(sm3);
        }
        else
        {
            float x, x2, mn, sm2, sm3;
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, ++X, ++Y)
                {
                    mn = sm2 = sm3 = 0.0f;
                    for (size_t l=0u; l<L-1; ++l, X+=K) { mn += *X; }
                    mn += *X; mn *= den;
                    for (size_t l=0u; l<L-1; ++l, X-=K) { x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; }
                    x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2;
                    *Y = w * sm3 / (sm2*sqrtf(sm2));
                }
            }
        }
    }

    return 0;
}


int skewness_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in skewness_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / L;
    const double w = (biased) ? sqrt(L) : L*sqrt(L-1)/(L-2);
    
    if (N==0u) {}
    else if (L<3) { fprintf(stderr,"error in skewness_d: L must be > 2\n"); return 1; }
    else if (L==N)
    {
        double x, x2, mn = 0.0, sm2 = 0.0, sm3 = 0.0;
        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=0u; l<L; ++l) { x = *--X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; }
        *Y = w * sm3 / (sm2*sqrt(sm2));
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            double x, x2, mn, sm2, sm3;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mn = sm2 = sm3 = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                X -= L;
                mn *= den;
                for (size_t l=0u; l<L; ++l, ++X) { x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; }
                *Y = w * sm3 / (sm2*sqrt(sm2));
            }
        }
        else if (G==1)
        {
            double x, x2, *sm2, *sm3;
            if (!(sm2=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm3=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += *X; }
            }
            X -= N;
            for (size_t v=0u; v<V; ++v, ++Y) { *Y *= den; }
            Y -= V;
            for (size_t l=0u; l<L; ++l, sm2-=V, sm3-=V, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++sm2, ++sm3, ++Y) { x = *X - *Y; x2 = x*x; *sm2 += x2; *sm3 += x*x2; }
            }
            for (size_t v=0u; v<V; ++v, ++sm2, ++sm3, ++Y) { *Y = w * *sm3 / (*sm2*sqrt(*sm2)); }
            sm2 -= V; sm3 -= V;
            free(sm2); free(sm3);
        }
        else
        {
            double x, x2, mn, sm2, sm3;
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, ++X, ++Y)
                {
                    mn = sm2 = sm3 = 0.0;
                    for (size_t l=0u; l<L-1; ++l, X+=K) { mn += *X; }
                    mn += *X; mn *= den;
                    for (size_t l=0u; l<L-1; ++l, X-=K) { x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; }
                    x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2;
                    *Y = w * sm3 / (sm2*sqrt(sm2));
                }
            }
        }
    }

    return 0;
}


int skewness_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in skewness_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / L;
    const float w = (biased) ? sqrtf(L) : L*sqrtf(L-1)/(L-2);
    float xr, xi, x2r, x2i, xrr, xii, xri, den3;

    if (N==0u) {}
    else if (L<3) { fprintf(stderr,"error in skewness_c: L must be > 2\n"); return 1; }
    else if (L==N)
    {
        float mnr = 0.0f, mni = 0.0f, sm2 = 0.0f, sm3r = 0.0f, sm3i = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
        mnr *= den; mni *= den;
        X -= 2*L;
        for (size_t l=0u; l<L; ++l, ++X)
        {
            xr = *X - mnr; xi = *++X - mni;
            xrr = xr*xr; xii = xi*xi; xri = xr*xi;
            x2r = xrr - xii; x2i = xri + xri;
            sm2 += xrr + xii;
            sm3r += xr*x2r - xi*x2i; sm3i += xr*x2i + xi*x2r;
        }
        den3 = w / (sm2*sqrtf(sm2));
        *Y = sm3r * den3; *++Y = sm3i * den3;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            float mnr, mni, sm2, sm3r, sm3i;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mnr = mni = sm2 = sm3r = sm3i = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
                mnr *= den; mni *= den;
                X -= 2*L;
                for (size_t l=0u; l<L; ++l, ++X)
                {
                    xr = *X - mnr; xi = *++X - mni;
                    xrr = xr*xr; xii = xi*xi; xri = xr*xi;
                    x2r = xrr - xii; x2i = xri + xri;
                    sm2 += xrr + xii;
                    sm3r += xr*x2r - xi*x2i; sm3i += xr*x2i + xi*x2r;
                }
                den3 = w / (sm2*sqrtf(sm2));
                *Y = sm3r * den3; *++Y = sm3i * den3;
            }
        }
        else
        {
            float mnr, mni, sm2, sm3r, sm3i;
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2, ++Y)
                {
                    mnr = mni = sm2 = sm3r = sm3i = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=2*K-1) { mnr += *X; mni += *++X; }
                    mnr *= den; mni *= den;
                    X -= 2*K*L;
                    for (size_t l=0u; l<L; ++l, X+=2*K-1)
                    {
                        xr = *X - mnr; xi = *++X - mni;
                        xrr = xr*xr; xii = xi*xi; xri = xr*xi;
                        x2r = xrr - xii; x2i = xri + xri;
                        sm2 += xrr + xii;
                        sm3r += xr*x2r - xi*x2i; sm3i += xr*x2i + xi*x2r;
                    }
                    den3 = w / (sm2*sqrtf(sm2));
                    *Y = sm3r * den3; *++Y = sm3i * den3;
                }
            }
        }
    }
    
    return 0;
}


int skewness_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in skewness_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / L;
    const double w = (biased) ? sqrt(L) : L*sqrt(L-1)/(L-2);
    double xr, xi, x2r, x2i, xrr, xii, xri, den3;

    if (N==0u) {}
    else if (L<3) { fprintf(stderr,"error in skewness_z: L must be > 2\n"); return 1; }
    else if (L==N)
    {
        double mnr = 0.0, mni = 0.0, sm2 = 0.0, sm3r = 0.0, sm3i = 0.0;
        for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
        mnr *= den; mni *= den;
        X -= 2*L;
        for (size_t l=0u; l<L; ++l, ++X)
        {
            xr = *X - mnr; xi = *++X - mni;
            xrr = xr*xr; xii = xi*xi; xri = xr*xi;
            x2r = xrr - xii; x2i = xri + xri;
            sm2 += xrr + xii;
            sm3r += xr*x2r - xi*x2i; sm3i += xr*x2i + xi*x2r;
        }
        den3 = w / (sm2*sqrt(sm2));
        *Y = sm3r * den3; *++Y = sm3i * den3;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            double mnr, mni, sm2, sm3r, sm3i;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mnr = mni = sm2 = sm3r = sm3i = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
                mnr *= den; mni *= den;
                X -= 2*L;
                for (size_t l=0u; l<L; ++l, ++X)
                {
                    xr = *X - mnr; xi = *++X - mni;
                    xrr = xr*xr; xii = xi*xi; xri = xr*xi;
                    x2r = xrr - xii; x2i = xri + xri;
                    sm2 += xrr + xii;
                    sm3r += xr*x2r - xi*x2i; sm3i += xr*x2i + xi*x2r;
                }
                den3 = w / (sm2*sqrt(sm2));
                *Y = sm3r * den3; *++Y = sm3i * den3;
            }
        }
        else
        {
            double mnr, mni, sm2, sm3r, sm3i;
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2, ++Y)
                {
                    mnr = mni = sm2 = sm3r = sm3i = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=2*K-1) { mnr += *X; mni += *++X; }
                    mnr *= den; mni *= den;
                    X -= 2*K*L;
                    for (size_t l=0u; l<L; ++l, X+=2*K-1)
                    {
                        xr = *X - mnr; xi = *++X - mni;
                        xrr = xr*xr; xii = xi*xi; xri = xr*xi;
                        x2r = xrr - xii; x2i = xri + xri;
                        sm2 += xrr + xii;
                        sm3r += xr*x2r - xi*x2i; sm3i += xr*x2i + xi*x2r;
                    }
                    den3 = w / (sm2*sqrt(sm2));
                    *Y = sm3r * den3; *++Y = sm3i * den3;
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
