//Vec2scalar (reduction) operation.
//Gets the kurtosis for each vector in X along dim.

//For complex case, output is complex.
//I follow the Octave convention for complex kurtosis (but see literature for other ideas later).

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int kurtosis_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int kurtosis_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int kurtosis_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int kurtosis_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);


int kurtosis_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in kurtosis_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / L;

    if (N==0) {}
    else if (L<4) { fprintf(stderr,"error in kurtosis_s: L must be > 3\n"); return 1; }
    else if (L==N)
    {
        float x, x2, sm2 = 0.0f, sm4 = 0.0f;
        if (L<7000)
        {
            *Y = 0.0f;
            for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
            *Y *= den;
            for (size_t l=0; l<L; ++l) { x = *--X - *Y; x2 = x*x; sm2 += x2; sm4 += x2*x2; }
        }
        else
        {
            const float o = 1.0f;
            *Y = cblas_sdot((int)L,X,1,&o,0) * den;
            for (size_t l=0; l<L; ++l, ++X) { x = *X - *Y; x2 = x*x; sm2 += x2; sm4 += x2*x2; }
        }
        *Y = L * sm4 / (sm2*sm2);
        if (!biased) { *Y =  3.0f + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            float x, x2, sm2, sm4;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                *Y = sm2 = sm4 = 0.0f;
                for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
                *Y *= den; X -= L;
                for (size_t l=0; l<L; ++l, ++X) { x = *X - *Y; x2 = x*x; sm2 += x2; sm4 += x2*x2; }
                *Y = L * sm4 / (sm2*sm2);
                if (!biased) { *Y =  3.0f + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
            }
        }
        else if (G==1)
        {
            const float z = 0.0f;
            float x, x2, *sm2, *sm4;
            if (!(sm2=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm4=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with calloc. "); perror("calloc"); return 1; }
            cblas_scopy((int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += *X; }
            }
            X -= N;
            cblas_sscal((int)V,den,Y,1);
            for (size_t l=0; l<L; ++l, sm2-=V, sm4-=V, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++sm2, ++sm4, ++Y) { x = *X - *Y; x2 = x*x; *sm2 += x2; *sm4 += x2*x2; }
            }
            for (size_t v=0; v<V; ++v, ++sm2, ++sm4, ++Y)
            {
                *Y = L * *sm4 / (*sm2**sm2);
                if (!biased) { *Y =  3.0f + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
            }
            sm2 -= V; sm4 -= V; free(sm2); free(sm4);
        }
        else
        {
            float x, x2, sm2, sm4, *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in kurtosis_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_scopy((int)L,X,(int)K,X1,1);
                    *Y = sm2 = sm4 = 0.0f;
                    for (size_t l=0; l<L; ++l, ++X1) { *Y += *X1; }
                    *Y *= den;
                    for (size_t l=0; l<L; ++l) { x = *X1-- - *Y; x2 = x*x; sm2 += x2; sm4 += x2*x2; }
                    *Y = L * sm4 / (sm2*sm2);
                    if (!biased) { *Y =  3.0f + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int kurtosis_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in kurtosis_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / L;

    if (N==0) {}
    else if (L<4) { fprintf(stderr,"error in kurtosis_d: L must be > 3\n"); return 1; }
    else if (L==N)
    {
        double x, x2, sm2 = 0.0, sm4 = 0.0;
        if (L<7000)
        {
            *Y = 0.0;
            for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
            *Y *= den;
            for (size_t l=0; l<L; ++l) { x = *--X - *Y; x2 = x*x; sm2 += x2; sm4 += x2*x2; }
        }
        else
        {
            const double o = 1.0;
            *Y = cblas_ddot((int)L,X,1,&o,0) * den;
            for (size_t l=0; l<L; ++l, ++X) { x = *X - *Y; x2 = x*x; sm2 += x2; sm4 += x2*x2; }
        }
        *Y = L * sm4 / (sm2*sm2);
        if (!biased) { *Y =  3.0 + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            double x, x2, sm2, sm4;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                *Y = sm2 = sm4 = 0.0;
                for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
                *Y *= den; X -= L;
                for (size_t l=0; l<L; ++l, ++X) { x = *X - *Y; x2 = x*x; sm2 += x2; sm4 += x2*x2; }
                *Y = L * sm4 / (sm2*sm2);
                if (!biased) { *Y =  3.0 + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
            }
        }
        else if (G==1)
        {
            const double z = 0.0;
            double x, x2, *sm2, *sm4;
            if (!(sm2=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with calloc. "); perror("calloc"); return 1; }
            if (!(sm4=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with calloc. "); perror("calloc"); return 1; }
            cblas_dcopy((int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += *X; }
            }
            X -= N;
            cblas_dscal((int)V,den,Y,1);
            for (size_t l=0; l<L; ++l, sm2-=V, sm4-=V, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++sm2, ++sm4, ++Y) { x = *X - *Y; x2 = x*x; *sm2 += x2; *sm4 += x2*x2; }
            }
            for (size_t v=0; v<V; ++v, ++sm2, ++sm4, ++Y)
            {
                *Y = L * *sm4 / (*sm2**sm2);
                if (!biased) { *Y =  3.0 + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
            }
            sm2 -= V; sm4 -= V; free(sm2); free(sm4);
        }
        else
        {
            double x, x2, sm2, sm4, *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in kurtosis_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_dcopy((int)L,X,(int)K,X1,1);
                    *Y = sm2 = sm4 = 0.0;
                    for (size_t l=0; l<L; ++l, ++X1) { *Y += *X1; }
                    *Y *= den;
                    for (size_t l=0; l<L; ++l) { x = *X1-- - *Y; x2 = x*x; sm2 += x2; sm4 += x2*x2; }
                    *Y = L * sm4 / (sm2*sm2);
                    if (!biased) { *Y =  3.0 + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3)); }
                }
            }
            free(X1);
        }
    }

    return 0;
}


int kurtosis_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in kurtosis_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / L;
    float xr, xi, x2r, x2i, x3r, x3i, xrr, xii, xri, den4;

    if (N==0) {}
    else if (L<4) { fprintf(stderr,"error in kurtosis_c: L must be > 3\n"); return 1; }
    else if (L==N)
    {
        float mnr = 0.0f, mni = 0.0f, sm2 = 0.0f, sm4r = 0.0f, sm4i = 0.0f;
        for (size_t l=0; l<L; ++l) { mnr += *X++; mni += *X++; }
        mnr *= den; mni *= den;
        X -= 2*L;
        for (size_t l=0; l<L; ++l)
        {
            xr = *X++ - mnr; xi = *X++ - mni;
            xrr = xr*xr; xii = xi*xi; xri = xr*xi;
            x2r = xrr - xii; x2i = xri + xri;
            sm2 += xrr + xii;
            x3r = xr*x2r - xi*x2i; x3i = xr*x2i + xi*x2r;
            sm4r += xr*x3r - xi*x3i; sm4i += xr*x3i + xi*x3r;
        }
        den4 = L / (sm2*sm2);
        *Y++ = sm4r * den4; *Y = sm4i * den4;
        if (!biased)
        {
            *Y-- *= (L+1)*(L-1) / (float)((L-2)*(L-3));
            *Y = 3.0f + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            float mnr, mni, sm2, sm4r, sm4i;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                mnr = mni = sm2 = sm4r = sm4i = 0.0f;
                for (size_t l=0; l<L; ++l) { mnr += *X++; mni += *X++; }
                mnr *= den; mni *= den;
                X -= 2*L;
                for (size_t l=0; l<L; ++l)
                {
                    xr = *X++ - mnr; xi = *X++ - mni;
                    xrr = xr*xr; xii = xi*xi; xri = xr*xi;
                    x2r = xrr - xii; x2i = xri + xri;
                    sm2 += xrr + xii;
                    x3r = xr*x2r - xi*x2i; x3i = xr*x2i + xi*x2r;
                    sm4r += xr*x3r - xi*x3i; sm4i += xr*x3i + xi*x3r;
                }
                den4 = L / (sm2*sm2);
                *Y++ = sm4r * den4; *Y = sm4i * den4;
                if (!biased)
                {
                    --Y;
                    *Y = 3.0f + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
                    *++Y *= (L+1)*(L-1) / (float)((L-2)*(L-3));
                }
            }
        }
        else
        {
            float mnr, mni, sm2, sm4r, sm4i, *X1;
            if (!(X1=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in kurtosis_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X+=2, X1-=2*L, ++Y)
                {
                    cblas_ccopy((int)L,X,(int)K,X1,1);
                    mnr = mni = sm2 = sm4r = sm4i = 0.0f;
                    for (size_t l=0; l<L; ++l) { mnr += *X1++; mni += *X1++; }
                    mnr *= den; mni *= den;
                    X1 -= 2*L;
                    for (size_t l=0; l<L; ++l)
                    {
                        xr = *X1++ - mnr; xi = *X1++ - mni;
                        xrr = xr*xr; xii = xi*xi; xri = xr*xi;
                        x2r = xrr - xii; x2i = xri + xri;
                        sm2 += xrr + xii;
                        x3r = xr*x2r - xi*x2i; x3i = xr*x2i + xi*x2r;
                        sm4r += xr*x3r - xi*x3i; sm4i += xr*x3i + xi*x3r;
                    }
                    den4 = L / (sm2*sm2);
                    *Y++ = sm4r * den4; *Y = sm4i * den4;
                    if (!biased)
                    {
                        --Y;
                        *Y = 3.0f + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
                        *++Y *= (L+1)*(L-1) / (float)((L-2)*(L-3));
                    }
                }
            }
            free(X1);
        }
    }
    
    return 0;
}


int kurtosis_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in kurtosis_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / L;
    double xr, xi, x2r, x2i, x3r, x3i, xrr, xii, xri, den4;

    if (N==0) {}
    else if (L<4) { fprintf(stderr,"error in kurtosis_z: L must be > 3\n"); return 1; }
    else if (L==N)
    {
        double mnr = 0.0, mni = 0.0, sm2 = 0.0, sm4r = 0.0, sm4i = 0.0;
        for (size_t l=0; l<L; ++l) { mnr += *X++; mni += *X++; }
        mnr *= den; mni *= den;
        X -= 2*L;
        for (size_t l=0; l<L; ++l)
        {
            xr = *X++ - mnr; xi = *X++ - mni;
            xrr = xr*xr; xii = xi*xi; xri = xr*xi;
            x2r = xrr - xii; x2i = xri + xri;
            sm2 += xrr + xii;
            x3r = xr*x2r - xi*x2i; x3i = xr*x2i + xi*x2r;
            sm4r += xr*x3r - xi*x3i; sm4i += xr*x3i + xi*x3r;
        }
        den4 = L / (sm2*sm2);
        *Y++ = sm4r * den4; *Y = sm4i * den4;
        if (!biased)
        {
            *Y-- *= (L+1)*(L-1) / (double)((L-2)*(L-3));
            *Y = 3.0 + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            double mnr, mni, sm2, sm4r, sm4i;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                mnr = mni = sm2 = sm4r = sm4i = 0.0;
                for (size_t l=0; l<L; ++l) { mnr += *X++; mni += *X++; }
                mnr *= den; mni *= den;
                X -= 2*L;
                for (size_t l=0; l<L; ++l)
                {
                    xr = *X++ - mnr; xi = *X++ - mni;
                    xrr = xr*xr; xii = xi*xi; xri = xr*xi;
                    x2r = xrr - xii; x2i = xri + xri;
                    sm2 += xrr + xii;
                    x3r = xr*x2r - xi*x2i; x3i = xr*x2i + xi*x2r;
                    sm4r += xr*x3r - xi*x3i; sm4i += xr*x3i + xi*x3r;
                }
                den4 = L / (sm2*sm2);
                *Y++ = sm4r * den4; *Y = sm4i * den4;
                if (!biased)
                {
                    --Y;
                    *Y = 3.0 + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
                    *++Y *= (L+1)*(L-1) / (double)((L-2)*(L-3));
                }
            }
        }
        else
        {
            double mnr, mni, sm2, sm4r, sm4i, *X1;
            if (!(X1=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in kurtosis_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X+=2, X1-=2*L, ++Y)
                {
                    cblas_zcopy((int)L,X,(int)K,X1,1);
                    mnr = mni = sm2 = sm4r = sm4i = 0.0;
                    for (size_t l=0; l<L; ++l) { mnr += *X1++; mni += *X1++; }
                    mnr *= den; mni *= den;
                    X1 -= 2*L;
                    for (size_t l=0; l<L; ++l)
                    {
                        xr = *X1++ - mnr; xi = *X1++ - mni;
                        xrr = xr*xr; xii = xi*xi; xri = xr*xi;
                        x2r = xrr - xii; x2i = xri + xri;
                        sm2 += xrr + xii;
                        x3r = xr*x2r - xi*x2i; x3i = xr*x2i + xi*x2r;
                        sm4r += xr*x3r - xi*x3i; sm4i += xr*x3i + xi*x3r;
                    }
                    den4 = L / (sm2*sm2);
                    *Y++ = sm4r * den4; *Y = sm4i * den4;
                    if (!biased)
                    {
                        --Y;
                        *Y = 3.0 + (*Y*(L+1)-3*(L-1)) * (L-1)/((L-2)*(L-3));
                        *++Y *= (L+1)*(L-1) / (double)((L-2)*(L-3));
                    }
                }
            }
            free(X1);
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
