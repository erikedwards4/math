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

int skewness_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
int skewness_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
int skewness_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
int skewness_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);


int skewness_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased)
{
    if (dim>3u) { fprintf(stderr,"error in skewness_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = 1.0f / (float)L;
    const float w = (biased) ? sqrtf((float)L) : (float)L*sqrtf((float)(L-1u))/(float)(L-2u);

    if (N==0u) {}
    else if (L<3u) { fprintf(stderr,"error in skewness_s: L must be > 2\n"); return 1; }
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
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float x, x2, mn, sm2, sm3;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                mn = sm2 = sm3 = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                X -= L;
                mn *= den;
                for (size_t l=0u; l<L; ++l, ++X) { x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; }
                *Y = w * sm3 / (sm2*sqrtf(sm2));
            }
        }
        else if (G==1u)
        {
            float x, x2, *sm2, *sm3;
            if (!(sm2=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm3=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in skewness_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += *X; }
            }
            X -= N;
            for (size_t v=V; v>0u; --v, ++Y) { *Y *= den; }
            Y -= V;
            for (size_t l=0u; l<L; ++l, sm2-=V, sm3-=V, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++sm2, ++sm3, ++Y) { x = *X - *Y; x2 = x*x; *sm2 += x2; *sm3 += x*x2; }
            }
            for (size_t v=V; v>0u; --v, ++sm2, ++sm3, ++Y) { *Y = w * *sm3 / (*sm2*sqrtf(*sm2)); }
            sm2 -= V; sm3 -= V;
            free(sm2); free(sm3);
        }
        else
        {
            float x, x2, mn, sm2, sm3;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X, ++Y)
                {
                    mn = sm2 = sm3 = 0.0f;
                    for (size_t l=0u; l<L-1u; ++l, X+=K) { mn += *X; }
                    mn += *X; mn *= den;
                    for (size_t l=0u; l<L-1u; ++l, X-=K) { x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; }
                    x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2;
                    *Y = w * sm3 / (sm2*sqrtf(sm2));
                }
            }
        }
    }

    return 0;
}


int skewness_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased)
{
    if (dim>3u) { fprintf(stderr,"error in skewness_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = 1.0 / (double)L;
    const double w = (biased) ? sqrt((double)L) : (double)L*sqrt((double)(L-1u))/(double)(L-2u);
    
    if (N==0u) {}
    else if (L<3u) { fprintf(stderr,"error in skewness_d: L must be > 2\n"); return 1; }
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
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double x, x2, mn, sm2, sm3;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                mn = sm2 = sm3 = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
                X -= L;
                mn *= den;
                for (size_t l=0u; l<L; ++l, ++X) { x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; }
                *Y = w * sm3 / (sm2*sqrt(sm2));
            }
        }
        else if (G==1u)
        {
            double x, x2, *sm2, *sm3;
            if (!(sm2=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm3=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in skewness_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += *X; }
            }
            X -= N;
            for (size_t v=V; v>0u; --v, ++Y) { *Y *= den; }
            Y -= V;
            for (size_t l=0u; l<L; ++l, sm2-=V, sm3-=V, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++sm2, ++sm3, ++Y) { x = *X - *Y; x2 = x*x; *sm2 += x2; *sm3 += x*x2; }
            }
            for (size_t v=V; v>0u; --v, ++sm2, ++sm3, ++Y) { *Y = w * *sm3 / (*sm2*sqrt(*sm2)); }
            sm2 -= V; sm3 -= V;
            free(sm2); free(sm3);
        }
        else
        {
            double x, x2, mn, sm2, sm3;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X, ++Y)
                {
                    mn = sm2 = sm3 = 0.0;
                    for (size_t l=0u; l<L-1u; ++l, X+=K) { mn += *X; }
                    mn += *X; mn *= den;
                    for (size_t l=0u; l<L-1u; ++l, X-=K) { x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; }
                    x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2;
                    *Y = w * sm3 / (sm2*sqrt(sm2));
                }
            }
        }
    }

    return 0;
}


int skewness_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased)
{
    if (dim>3u) { fprintf(stderr,"error in skewness_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = 1.0f / (float)L;
    const float w = (biased) ? sqrtf((float)L) : (float)L*sqrtf((float)(L-1u))/(float)(L-2u);
    float xr, xi, x2r, x2i, xrr, xii, xri, den3;

    if (N==0u) {}
    else if (L<3u) { fprintf(stderr,"error in skewness_c: L must be > 2\n"); return 1; }
    else if (L==N)
    {
        float mnr = 0.0f, mni = 0.0f, sm2 = 0.0f, sm3r = 0.0f, sm3i = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
        mnr *= den; mni *= den;
        X -= 2u*L;
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
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float mnr, mni, sm2, sm3r, sm3i;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                mnr = mni = sm2 = sm3r = sm3i = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
                mnr *= den; mni *= den;
                X -= 2u*L;
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
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    mnr = mni = sm2 = sm3r = sm3i = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u) { mnr += *X; mni += *++X; }
                    mnr *= den; mni *= den;
                    X -= 2u*K*L;
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u)
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


int skewness_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased)
{
    if (dim>3u) { fprintf(stderr,"error in skewness_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = 1.0 / (double)L;
    const double w = (biased) ? sqrt((double)L) : (double)L*sqrt((double)(L-1u))/(double)(L-2u);
    double xr, xi, x2r, x2i, xrr, xii, xri, den3;

    if (N==0u) {}
    else if (L<3u) { fprintf(stderr,"error in skewness_z: L must be > 2\n"); return 1; }
    else if (L==N)
    {
        double mnr = 0.0, mni = 0.0, sm2 = 0.0, sm3r = 0.0, sm3i = 0.0;
        for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
        mnr *= den; mni *= den;
        X -= 2u*L;
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
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double mnr, mni, sm2, sm3r, sm3i;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                mnr = mni = sm2 = sm3r = sm3i = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { mnr += *X; mni += *++X; }
                mnr *= den; mni *= den;
                X -= 2u*L;
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
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    mnr = mni = sm2 = sm3r = sm3i = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u) { mnr += *X; mni += *++X; }
                    mnr *= den; mni *= den;
                    X -= 2u*K*L;
                    for (size_t l=0u; l<L; ++l, X+=2u*K-1u)
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
