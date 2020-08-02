//Vec2scalar (reduction) operation.
//Gets variance for each vector in X along dim.
//For complex case, output is real.

//This originally worked in place for some conditions, but those are deleted now.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int var_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int var_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int var_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int var_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);


int var_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in var_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f/L, den2 = (biased) ? den : 1.0f/(L-1);

    if (N==0) {}
    else if (L==1)
    {
        const float z = 0.0f;
        cblas_scopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        float x, sm2 = 0.0f;
        if (L<7000)
        {
            *Y = 0.0f;
            for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
            *Y *= den;
            for (size_t l=0; l<L; ++l) { x = *--X - *Y; sm2 += x*x; }
        }
        else
        {
            const float o = 1.0f;
            *Y = cblas_sdot((int)L,X,1,&o,0) * den;
            for (size_t l=0; l<L; ++l, ++X) { x = *X - *Y; sm2 += x*x; }
        }
        *Y = sm2 * den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            float x, sm2;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                *Y = sm2 = 0.0f;
                for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
                *Y *= den; X -= L;
                for (size_t l=0; l<L; ++l, ++X) { x = *X - *Y; sm2 += x*x; }
                *Y = sm2 * den2;
            }
        }
        else if (G==1)
        {
            const float z = 0.0f;
            float x, *mn;
            if (!(mn=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in var_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0; l<L; ++l, mn-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++mn) { *mn += *X; }
            }
            X -= N;
            cblas_sscal((int)V,den,mn,1);
            cblas_scopy((int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, mn-=V, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++mn, ++Y) { x = *X - *mn; *Y += x*x; }
            }
            cblas_sscal((int)V,den2,Y,1);
            free(mn);
        }
        else
        {
            float x, sm2, *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in var_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_scopy((int)L,X,(int)K,X1,1);
                    *Y = sm2 = 0.0f;
                    for (size_t l=0; l<L; ++l, ++X1) { *Y += *X1; }
                    *Y *= den;
                    for (size_t l=0; l<L; ++l) { x = *--X1 - *Y; sm2 += x*x; }
                    *Y = sm2 * den2;
                }
            }
        }
    }

    return 0;
}


int var_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in var_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0/L, den2 = (biased) ? den : 1.0/(L-1);
    
    if (N==0) {}
    else if (L==1)
    {
        const double z = 0.0;
        cblas_dcopy((int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        double x, sm2 = 0.0;
        if (L<7000)
        {
            *Y = 0.0;
            for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
            *Y *= den;
            for (size_t l=0; l<L; ++l) { x = *--X - *Y; sm2 += x*x; }
        }
        else
        {
            const double o = 1.0;
            *Y = cblas_ddot((int)L,X,1,&o,0) * den;
            for (size_t l=0; l<L; ++l, ++X) { x = *X - *Y; sm2 += x*x; }
        }
        *Y = sm2 * den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            double x, sm2;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                *Y = sm2 = 0.0;
                for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
                *Y *= den; X -= L;
                for (size_t l=0; l<L; ++l, ++X) { x = *X - *Y; sm2 += x*x; }
                *Y = sm2 * den2;
            }
        }
        else if (G==1)
        {
            const double z = 0.0;
            double x, *mn;
            if (!(mn=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in var_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0; l<L; ++l, mn-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++mn) { *mn += *X; }
            }
            X -= N;
            cblas_dscal((int)V,den,mn,1);
            cblas_dcopy((int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, mn-=V, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++mn, ++Y) { x = *X - *mn; *Y += x*x; }
            }
            cblas_dscal((int)V,den2,Y,1);
            free(mn);
        }
        else
        {
            double x, sm2, *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in var_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_dcopy((int)L,X,(int)K,X1,1);
                    *Y = sm2 = 0.0;
                    for (size_t l=0; l<L; ++l, ++X1) { *Y += *X1; }
                    *Y *= den;
                    for (size_t l=0; l<L; ++l) { x = *--X1 - *Y; sm2 += x*x; }
                    *Y = sm2 * den2;
                }
            }
        }
    }

    return 0;
}


int var_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in var_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f/L, den2 = (biased) ? den : 1.0f/(L-1);
    float xr, xi;

    if (N==0) {}
    else if (L==1)
    {
        const float z = 0.0f;
        cblas_scopy(2*(int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        float mnr = 0.0f, mni = 0.0f, sm2 = 0.0f;
        for (size_t l=0; l<L; ++l, ++X) { mnr += *X++; mni += *X; }
        mnr *= den; mni *= den;
        X -= 2*L;
        for (size_t l=0; l<L; ++l, ++X) { xr = *X++ - mnr; xi = *X - mni; sm2 += xr*xr + xi*xi; }
        *Y = sm2 * den2;
        // else //this is faster for L>12000, but not in-place and not usable for skewness, etc.
        // {
        //     const float den[2] = {-1.0f/L,0.0f};
        //     cblas_caxpy((int)L,mn,den,0,X,1);
        //     cblas_cdotc_sub((int)L,X,1,X,1,(_Complex float *)mn);
        //     *Y = mn[0] * den2;
        // }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            float mnr, mni, sm2;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                mnr = mni = sm2 = 0.0f;
                for (size_t l=0; l<L; ++l, ++X) { mnr += *X++; mni += *X; }
                mnr *= den; mni *= den;
                X -= 2*L;
                for (size_t l=0; l<L; ++l, ++X) { xr = *X++ - mnr; xi = *X - mni; sm2 += xr*xr + xi*xi; }
                *Y = sm2 * den2;
            }
        }
        else if (G==1)
        {
            const float z = 0.0f;
            float *mn;
            if (!(mn=(float *)calloc(2*V,sizeof(float)))) { fprintf(stderr,"error in var_c: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0; l<L; ++l, mn-=2*V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++mn) { *mn++ += *X++; *mn += *X; }
            }
            X -= 2*N;
            cblas_sscal(2*(int)V,den,mn,1);
            cblas_scopy(2*(int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, mn-=2*V, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++mn, ++Y) { xr = *X++ - *mn++; xi = *X - *mn; *Y += xr*xr + xi*xi; }
            }
            cblas_sscal((int)V,den2,Y,1);
            free(mn);
        }
        else
        {
            float mnr, mni, sm2, *X1;
            if (!(X1=(float *)malloc(2*L*sizeof(float)))) { fprintf(stderr,"error in var_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X+=2, X1-=2*L, ++Y)
                {
                    cblas_ccopy((int)L,X,(int)K,X1,1);
                    mnr = mni = sm2 = 0.0f;
                    for (size_t l=0; l<L; l++, ++X1) { mnr += *X1++; mni += *X1; }
                    mnr *= den; mni *= den;
                    X1 -= 2*L;
                    for (size_t l=0; l<L; ++l, ++X1) { xr = *X1++ - mnr; xi = *X1 - mni; *Y += xr*xr + xi*xi; }
                    *Y = sm2 * den2;
                }
            }
            free(X1);
        }
    }
    
    return 0;
}


int var_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in var_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0/L, den2 = (biased) ? den : 1.0/(L-1);
    double xr, xi;

    if (N==0) {}
    else if (L==1)
    {
        const double z = 0.0;
        cblas_dcopy(2*(int)N,&z,0,Y,1);
    }
    else if (L==N)
    {
        double mnr = 0.0, mni = 0.0, sm2 = 0.0;
        for (size_t l=0; l<L; ++l) { mnr += *X++; mni += *X++; }
        mnr *= den; mni *= den;
        X -= 2*L;
        for (size_t l=0; l<L; ++l) { xr = *X++ - mnr; xi = *X++ - mni; sm2 += xr*xr + xi*xi; }
        *Y = sm2 * den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            double mnr, mni, sm2;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                mnr = mni = sm2 = 0.0;
                for (size_t l=0; l<L; ++l) { mnr += *X++; mni += *X++; }
                mnr *= den; mni *= den;
                X -= 2*L;
                for (size_t l=0; l<L; ++l) { xr = *X++ - mnr; xi = *X++ - mni; sm2 += xr*xr + xi*xi; }
                *Y = sm2 * den2;
            }
        }
        else if (G==1)
        {
            const double z = 0.0;
            double *mn;
            if (!(mn=(double *)calloc(2*V,sizeof(double)))) { fprintf(stderr,"error in var_z: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0; l<L; ++l, mn-=2*V)
            {
                for (size_t v=0; v<V; ++v) { *mn++ += *X++; *mn++ += *X++; }
            }
            X -= 2*N;
            cblas_dscal(2*(int)V,den,mn,1);
            cblas_dcopy(2*(int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, mn-=2*V, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++Y) { xr = *X++ - *mn++; xi = *X++ - *mn++; *Y += xr*xr + xi*xi; }
            }
            cblas_dscal((int)V,den2,Y,1);
            free(mn);
        }
        else
        {
            double mnr, mni, sm2, *X1;
            if (!(X1=(double *)malloc(2*L*sizeof(double)))) { fprintf(stderr,"error in var_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X+=2, X1-=2*L, ++Y)
                {
                    cblas_zcopy((int)L,X,(int)K,X1,1);
                    mnr = mni = sm2 = 0.0;
                    for (size_t l=0; l<L; l++, ++X1) { mnr += *X1++; mni += *X1; }
                    mnr *= den; mni *= den;
                    X1 -= 2*L;
                    for (size_t l=0; l<L; ++l, ++X1) { xr = *X1++ - mnr; xi = *X1 - mni; *Y += xr*xr + xi*xi; }
                    *Y = sm2 * den2;
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
