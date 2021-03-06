//Vec2vec operation.
//Gets the first 4 statistical moments for each vector in X along dim.
//These are mean, var, skewness, kurtosis, with the same biased option as usual.

//For complex case, output is complex.
//I follow the Octave convention for complex moments (but see literature for other ideas later).

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int moments_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int moments_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int moments_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int moments_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);


int moments_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in moments_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Ly = 4;
    const float den = 1.0f / Lx, den2 = (biased) ? den : 1.0f/(Lx-1);
    const float w = (biased) ? sqrtf(Lx) : Lx*sqrtf(Lx-1)/(Lx-2);

    if (N==0u) {}
    else if (Lx<4) { fprintf(stderr,"error in moments_s: L must be > 3\n"); return 1; }
    else if (Lx==N)
    {
        float x, x2, mn = 0.0f, sm2 = 0.0f, sm3 = 0.0f, sm4 = 0.0f;
        mn = 0.0f;
        for (size_t l=0u; l<Lx; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=0u; l<Lx; ++l) { x = *--X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; sm4 += x2*x2; }
        *Y = mn;
        *++Y = sm2 * den2;
        *++Y = w * sm3 / (sm2*sqrtf(sm2));
        *++Y = Lx * sm4 / (sm2*sm2);
        if (!biased) { *Y =  3.0f + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3)); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            float x, x2, mn, sm2, sm3, sm4;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mn = sm2 = sm3 = sm4 = 0.0f;
                for (size_t l=0u; l<Lx; ++l, ++X) { mn += *X; }
                X -= Lx;
                mn *= den;
                for (size_t l=0u; l<Lx; ++l, ++X) { x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; sm4 += x2*x2; }
                *Y = mn;
                *++Y = sm2 * den2;
                *++Y = w * sm3 / (sm2*sqrtf(sm2));
                *++Y = Lx * sm4 / (sm2*sm2);
                if (!biased) { *Y =  3.0f + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3)); }
            }
        }
        else if (G==1)
        {
            float x, x2, *sm2, *sm3, *sm4;
            if (!(sm2=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in moments_s: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm3=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in moments_s: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm4=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in moments_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1u; l<Lx; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += *X; }
            }
            X -= N;
            for (size_t v=0u; v<V; ++v, ++Y) { *Y *= den; }
            Y -= V;
            for (size_t l=0u; l<Lx; ++l, sm2-=V, sm3-=V, sm4-=V, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++sm2, ++sm3, ++sm4, ++Y) { x = *X - *Y; x2 = x*x; *sm2 += x2; *sm3 += x*x2; *sm4 += x2*x2; }
            }
            Y += V;
            for (size_t v=0u; v<V; ++v, ++sm2, ++sm3, ++sm4, ++Y) { *Y = *sm2 * den2; }
            sm2 -= V; sm3 -= V; sm4 -= V;
            for (size_t v=0u; v<V; ++v, ++sm2, ++sm3, ++sm4, ++Y) { *Y = w * *sm3 / (*sm2*sqrtf(*sm2)); }
            sm2 -= V; sm3 -= V; sm4 -= V;
            for (size_t v=0u; v<V; ++v, ++sm2, ++sm3, ++sm4, ++Y)
            {
                *Y = Lx * *sm4 / (*sm2**sm2);
                if (!biased) { *Y =  3.0f + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3)); }
            }
            sm2 -= V; sm3 -= V; sm4 -= V;
            free(sm2); free(sm3); free(sm4);
        }
        else
        {
            float x, x2, mn, sm2, sm3, sm4;
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1), Y+=B*(Ly-1))
            {
                for (size_t b=0u; b<B; ++b, ++X, Y-=(Ly-1)*K-1)
                {
                    mn = sm2 = sm3 = sm4 = 0.0f;
                    for (size_t l=0u; l<Lx-1; ++l, X+=K) { mn += *X; }
                    mn += *X; mn *= den;
                    for (size_t l=0u; l<Lx-1; ++l, X-=K) { x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; sm4 += x2*x2; }
                    x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; sm4 += x2*x2;
                    *Y = mn; Y += K;
                    *Y = sm2 * den2; Y += K;
                    *Y = w * sm3 / (sm2*sqrtf(sm2)); Y += K;
                    *Y = Lx * sm4 / (sm2*sm2);
                    if (!biased) { *Y =  3.0f + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3)); }
                }
            }
        }
    }

    return 0;
}


int moments_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in moments_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Ly = 4;
    const double den = 1.0 / Lx, den2 = (biased) ? den : 1.0/(Lx-1);
    const double w = (biased) ? sqrt(Lx) : Lx*sqrt(Lx-1)/(Lx-2);

    if (N==0u) {}
    else if (Lx<4) { fprintf(stderr,"error in moments_d: L must be > 3\n"); return 1; }
    else if (Lx==N)
    {
        double x, x2, mn = 0.0, sm2 = 0.0, sm3 = 0.0, sm4 = 0.0;
        mn = 0.0;
        for (size_t l=0u; l<Lx; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=0u; l<Lx; ++l) { x = *--X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; sm4 += x2*x2; }
        *Y = mn;
        *++Y = sm2 * den2;
        *++Y = w * sm3 / (sm2*sqrt(sm2));
        *++Y = Lx * sm4 / (sm2*sm2);
        if (!biased) { *Y =  3.0 + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3)); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            double x, x2, mn, sm2, sm3, sm4;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mn = sm2 = sm3 = sm4 = 0.0;
                for (size_t l=0u; l<Lx; ++l, ++X) { mn += *X; }
                X -= Lx;
                mn *= den;
                for (size_t l=0u; l<Lx; ++l, ++X) { x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; sm4 += x2*x2; }
                *Y = mn;
                *++Y = sm2 * den2;
                *++Y = w * sm3 / (sm2*sqrt(sm2));
                *++Y = Lx * sm4 / (sm2*sm2);
                if (!biased) { *Y =  3.0 + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3)); }
            }
        }
        else if (G==1)
        {
            double x, x2, *sm2, *sm3, *sm4;
            if (!(sm2=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in moments_d: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm3=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in moments_d: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm4=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in moments_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1u; l<Lx; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += *X; }
            }
            X -= N;
            for (size_t v=0u; v<V; ++v, ++Y) { *Y *= den; }
            Y -= V;
            for (size_t l=0u; l<Lx; ++l, sm2-=V, sm3-=V, sm4-=V, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++sm2, ++sm3, ++sm4, ++Y) { x = *X - *Y; x2 = x*x; *sm2 += x2; *sm3 += x*x2; *sm4 += x2*x2; }
            }
            Y += V;
            for (size_t v=0u; v<V; ++v, ++sm2, ++sm3, ++sm4, ++Y) { *Y = *sm2 * den2; }
            sm2 -= V; sm3 -= V; sm4 -= V;
            for (size_t v=0u; v<V; ++v, ++sm2, ++sm3, ++sm4, ++Y) { *Y = w * *sm3 / (*sm2*sqrt(*sm2)); }
            sm2 -= V; sm3 -= V; sm4 -= V;
            for (size_t v=0u; v<V; ++v, ++sm2, ++sm3, ++sm4, ++Y)
            {
                *Y = Lx * *sm4 / (*sm2**sm2);
                if (!biased) { *Y =  3.0 + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3)); }
            }
            sm2 -= V; sm3 -= V; sm4 -= V;
            free(sm2); free(sm3); free(sm4);
        }
        else
        {
            double x, x2, mn, sm2, sm3, sm4;
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1), Y+=B*(Ly-1))
            {
                for (size_t b=0u; b<B; ++b, ++X, Y-=(Ly-1)*K-1)
                {
                    mn = sm2 = sm3 = sm4 = 0.0;
                    for (size_t l=0u; l<Lx-1; ++l, X+=K) { mn += *X; }
                    mn += *X; mn *= den;
                    for (size_t l=0u; l<Lx-1; ++l, X-=K) { x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; sm4 += x2*x2; }
                    x = *X - mn; x2 = x*x; sm2 += x2; sm3 += x*x2; sm4 += x2*x2;
                    *Y = mn; Y += K;
                    *Y = sm2 * den2; Y += K;
                    *Y = w * sm3 / (sm2*sqrt(sm2)); Y += K;
                    *Y = Lx * sm4 / (sm2*sm2);
                    if (!biased) { *Y =  3.0 + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3)); }
                }
            }
        }
    }

    return 0;
}


int moments_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in moments_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Ly = 4;
    const float den = 1.0f/Lx, den2 = (biased) ? den : 1.0f/(Lx-1);
    const float w = (biased) ? sqrtf(Lx) : Lx*sqrtf(Lx-1)/(Lx-2);
    float xr, xi, x2r, x2i, x3r, x3i, xrr, xii, xri, den3, den4;
    float mnr = 0.0f, mni = 0.0f, sm2 = 0.0f, sm3r = 0.0f, sm3i = 0.0f, sm4r = 0.0f, sm4i = 0.0f;

    if (N==0u) {}
    else if (Lx<4) { fprintf(stderr,"error in moments_c: Lx must be > 3\n"); return 1; }
    else if (Lx==N)
    {
        for (size_t l=0u; l<Lx; ++l, ++X) { mnr += *X; mni += *++X; }
        mnr *= den; mni *= den;
        X -= 2*Lx;
        for (size_t l=0u; l<Lx; ++l, ++X)
        {
            xr = *X - mnr; xi = *++X - mni;
            xrr = xr*xr; xii = xi*xi; xri = xr*xi;
            x2r = xrr - xii; x2i = xri + xri;
            sm2 += xrr + xii;
            x3r = xr*x2r - xi*x2i; x3i = xr*x2i + xi*x2r;
            sm3r += x3r; sm3i += x3i;
            sm4r += xr*x3r - xi*x3i; sm4i += xr*x3i + xi*x3r;
        }
        *Y = mnr; *++Y = mni;
        *++Y = sm2 * den2; *++Y = 0.0f;
        den3 = w / (sm2*sqrtf(sm2));
        *++Y = sm3r * den3; *++Y = sm3i * den3;
        den4 = Lx / (sm2*sm2);
        *++Y = sm4r * den4; *++Y = sm4i * den4;
        if (!biased)
        {
            *Y-- *= (Lx+1)*(Lx-1) / (float)((Lx-2)*(Lx-3));
            *Y = 3.0f + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3));
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mnr = mni = sm2 = sm3r = sm3i = sm4r = sm4i = 0.0f;
                for (size_t l=0u; l<Lx; ++l, ++X) { mnr += *X; mni += *++X; }
                mnr *= den; mni *= den;
                X -= 2*Lx;
                for (size_t l=0u; l<Lx; ++l, ++X)
                {
                    xr = *X - mnr; xi = *++X - mni;
                    xrr = xr*xr; xii = xi*xi; xri = xr*xi;
                    x2r = xrr - xii; x2i = xri + xri;
                    sm2 += xrr + xii;
                    x3r = xr*x2r - xi*x2i; x3i = xr*x2i + xi*x2r;
                    sm3r += x3r; sm3i += x3i;
                    sm4r += xr*x3r - xi*x3i; sm4i += xr*x3i + xi*x3r;
                }
                *Y = mnr; *++Y = mni;
                *++Y = sm2 * den2; *++Y = 0.0f;
                den3 = w / (sm2*sqrtf(sm2));
                *++Y = sm3r * den3; *++Y = sm3i * den3;
                den4 = Lx / (sm2*sm2);
                *++Y = sm4r * den4; *++Y = sm4i * den4;
                if (!biased)
                {
                    --Y;
                    *Y = 3.0f + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3));
                    *++Y *= (Lx+1)*(Lx-1) / (float)((Lx-2)*(Lx-3));
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(Lx-1), Y+=2*B*(Ly-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*Lx-2, Y-=2*(Ly-1)*K-1)
                {
                    mnr = mni = sm2 = sm3r = sm3i = sm4r = sm4i = 0.0f;
                    for (size_t l=0u; l<Lx; ++l, X+=2*K-1) { mnr += *X; mni += *++X; }
                    mnr *= den; mni *= den;
                    X -= 2*K*Lx;
                    for (size_t l=0u; l<Lx; ++l, X+=2*K-1)
                    {
                        xr = *X - mnr; xi = *++X - mni;
                        xrr = xr*xr; xii = xi*xi; xri = xr*xi;
                        x2r = xrr - xii; x2i = xri + xri;
                        sm2 += xrr + xii;
                        x3r = xr*x2r - xi*x2i; x3i = xr*x2i + xi*x2r;
                        sm3r += x3r; sm3i += x3i;
                        sm4r += xr*x3r - xi*x3i; sm4i += xr*x3i + xi*x3r;
                    }
                    *Y = mnr; *++Y = mni; Y += 2*K-1;
                    *Y = sm2 * den2; *++Y = 0.0f; Y += 2*K-1;
                    den3 = w / (sm2*sqrtf(sm2));
                    *Y = sm3r * den3; *++Y = sm3i * den3; Y += 2*K-1;
                    den4 = Lx / (sm2*sm2);
                    *Y = sm4r * den4; *++Y = sm4i * den4;
                    if (!biased)
                    {
                        --Y;
                        *Y = 3.0f + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3));
                        *++Y *= (Lx+1)*(Lx-1) / (float)((Lx-2)*(Lx-3));
                    }
                }
            }
        }
    }
    
    return 0;
}


int moments_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in moments_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Ly = 4;
    const double den = 1.0/Lx, den2 = (biased) ? den : 1.0/(Lx-1);
    const double w = (biased) ? sqrt(Lx) : Lx*sqrt(Lx-1)/(Lx-2);
    double xr, xi, x2r, x2i, x3r, x3i, xrr, xii, xri, den3, den4;
    double mnr = 0.0, mni = 0.0, sm2 = 0.0, sm3r = 0.0, sm3i = 0.0, sm4r = 0.0, sm4i = 0.0;

    if (N==0u) {}
    else if (Lx<4) { fprintf(stderr,"error in moments_z: Lx must be > 3\n"); return 1; }
    else if (Lx==N)
    {
        for (size_t l=0u; l<Lx; ++l, ++X) { mnr += *X; mni += *++X; }
        mnr *= den; mni *= den;
        X -= 2*Lx;
        for (size_t l=0u; l<Lx; ++l, ++X)
        {
            xr = *X - mnr; xi = *++X - mni;
            xrr = xr*xr; xii = xi*xi; xri = xr*xi;
            x2r = xrr - xii; x2i = xri + xri;
            sm2 += xrr + xii;
            x3r = xr*x2r - xi*x2i; x3i = xr*x2i + xi*x2r;
            sm3r += x3r; sm3i += x3i;
            sm4r += xr*x3r - xi*x3i; sm4i += xr*x3i + xi*x3r;
        }
        *Y = mnr; *++Y = mni;
        *++Y = sm2 * den2; *++Y = 0.0;
        den3 = w / (sm2*sqrt(sm2));
        *++Y = sm3r * den3; *++Y = sm3i * den3;
        den4 = Lx / (sm2*sm2);
        *++Y = sm4r * den4; *++Y = sm4i * den4;
        if (!biased)
        {
            *Y-- *= (Lx+1)*(Lx-1) / (double)((Lx-2)*(Lx-3));
            *Y = 3.0 + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3));
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mnr = mni = sm2 = sm3r = sm3i = sm4r = sm4i = 0.0;
                for (size_t l=0u; l<Lx; ++l, ++X) { mnr += *X; mni += *++X; }
                mnr *= den; mni *= den;
                X -= 2*Lx;
                for (size_t l=0u; l<Lx; ++l, ++X)
                {
                    xr = *X - mnr; xi = *++X - mni;
                    xrr = xr*xr; xii = xi*xi; xri = xr*xi;
                    x2r = xrr - xii; x2i = xri + xri;
                    sm2 += xrr + xii;
                    x3r = xr*x2r - xi*x2i; x3i = xr*x2i + xi*x2r;
                    sm3r += x3r; sm3i += x3i;
                    sm4r += xr*x3r - xi*x3i; sm4i += xr*x3i + xi*x3r;
                }
                *Y = mnr; *++Y = mni;
                *++Y = sm2 * den2; *++Y = 0.0;
                den3 = w / (sm2*sqrt(sm2));
                *++Y = sm3r * den3; *++Y = sm3i * den3;
                den4 = Lx / (sm2*sm2);
                *++Y = sm4r * den4; *++Y = sm4i * den4;
                if (!biased)
                {
                    --Y;
                    *Y = 3.0 + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3));
                    *++Y *= (Lx+1)*(Lx-1) / (double)((Lx-2)*(Lx-3));
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(Lx-1), Y+=2*B*(Ly-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*Lx-2, Y-=2*(Ly-1)*K-1)
                {
                    mnr = mni = sm2 = sm3r = sm3i = sm4r = sm4i = 0.0;
                    for (size_t l=0u; l<Lx; ++l, X+=2*K-1) { mnr += *X; mni += *++X; }
                    mnr *= den; mni *= den;
                    X -= 2*K*Lx;
                    for (size_t l=0u; l<Lx; ++l, X+=2*K-1)
                    {
                        xr = *X - mnr; xi = *++X - mni;
                        xrr = xr*xr; xii = xi*xi; xri = xr*xi;
                        x2r = xrr - xii; x2i = xri + xri;
                        sm2 += xrr + xii;
                        x3r = xr*x2r - xi*x2i; x3i = xr*x2i + xi*x2r;
                        sm3r += x3r; sm3i += x3i;
                        sm4r += xr*x3r - xi*x3i; sm4i += xr*x3i + xi*x3r;
                    }
                    *Y = mnr; *++Y = mni; Y += 2*K-1;
                    *Y = sm2 * den2; *++Y = 0.0; Y += 2*K-1;
                    den3 = w / (sm2*sqrt(sm2));
                    *Y = sm3r * den3; *++Y = sm3i * den3; Y += 2*K-1;
                    den4 = Lx / (sm2*sm2);
                    *Y = sm4r * den4; *++Y = sm4i * den4;
                    if (!biased)
                    {
                        --Y;
                        *Y = 3.0 + (*Y*(Lx+1)-3*(Lx-1)) * (Lx-1)/((Lx-2)*(Lx-3));
                        *++Y *= (Lx+1)*(Lx-1) / (double)((Lx-2)*(Lx-3));
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
