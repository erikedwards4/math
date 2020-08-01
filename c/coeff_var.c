//Vec2scalar (reduction) operation.
//Gets the coefficient of variation (std/mean) for each vector in X along dim.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int coeff_var_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);
int coeff_var_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased);


int coeff_var_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in coeff_var_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float ni = 1.0f/L, den = (biased) ? ni : 1.0f/(L-1);

    if (N==0) {}
    else if (L<2) { fprintf(stderr,"error in coeff_var_s: L must be > 1\n"); return 1; }
    else if (L==N)
    {
        float x, sm2 = 0.0f;
        if (L<7000)
        {
            *Y = 0.0f;
            for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
            *Y *= ni;
            for (size_t l=0; l<L; ++l) { x = *--X - *Y; sm2 += x*x; }
        }
        else
        {
            const float o = 1.0f;
            *Y = cblas_sdot((int)L,X,1,&o,0) * ni;
            for (size_t l=0; l<L; ++l, ++X) { x = *X - *Y; sm2 += x*x; }
        }
        *Y = sqrtf(sm2*den) / *Y;
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
                *Y *= ni; X -= L;
                for (size_t l=0; l<L; ++l, ++X) { x = *X - *Y; sm2 += x*x; }
                *Y = sqrtf(sm2*den) / *Y;
            }
        }
        else if (G==1)
        {
            const float z = 0.0f;
            float x, *mn;
            if (!(mn=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0; l<L; ++l, mn-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++mn) { *mn += *X; }
            }
            X -= N;
            cblas_sscal((int)V,ni,mn,1);
            cblas_scopy((int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, mn-=V, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++mn, ++Y) { x = *X - *mn; *Y += x*x; }
            }
            for (size_t v=0; v<V; ++v, ++mn, ++Y) { *Y = sqrtf(*Y*den) / *mn; }
            mn -= V; free(mn);
        }
        else
        {
            float x, sm2, *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_scopy((int)L,X,(int)K,X1,1);
                    *Y = sm2 = 0.0;
                    for (size_t l=0; l<L; ++l, ++X1) { *Y += *X1; }
                    *Y *= ni;
                    for (size_t l=0; l<L; ++l) { x = *--X1 - *Y; sm2 += x*x; }
                    *Y = sqrtf(sm2*den) / *Y;
                }
            }
        }
    }

    return 0;
}


int coeff_var_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in coeff_var_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double ni = 1.0/L, den = (biased) ? ni : 1.0/(L-1);

    if (N==0) {}
    else if (L<2) { fprintf(stderr,"error in coeff_var_d: L must be > 1\n"); return 1; }
    else if (L==N)
    {
        double x, sm2 = 0.0;
        if (L<7000)
        {
            *Y = 0.0;
            for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
            *Y *= ni;
            for (size_t l=0; l<L; ++l) { x = *--X - *Y; sm2 += x*x; }
        }
        else
        {
            const double o = 1.0;
            *Y = cblas_ddot((int)L,X,1,&o,0) * ni;
            for (size_t l=0; l<L; ++l, ++X) { x = *X - *Y; sm2 += x*x; }
        }
        *Y = sqrt(sm2*den);
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
                *Y *= ni; X -= L;
                for (size_t l=0; l<L; ++l, ++X) { x = *X - *Y; sm2 += x*x; }
                *Y = sqrt(sm2*den);
            }
        }
        else if (G==1)
        {
            const double z = 0.0;
            double x, *mn;
            if (!(mn=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0; l<L; ++l, mn-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++mn) { *mn += *X; }
            }
            X -= N;
            cblas_dscal((int)V,ni,mn,1);
            cblas_dcopy((int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, mn-=V, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++mn, ++Y) { x = *X - *mn; *Y += x*x; }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = sqrt(*Y*den); }
            free(mn);
        }
        else
        {
            double x, sm2, *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                    cblas_dcopy((int)L,X,(int)K,X1,1);
                    *Y = sm2 = 0.0;
                    for (size_t l=0; l<L; ++l, ++X1) { *Y += *X1; }
                    *Y *= ni;
                    for (size_t l=0; l<L; ++l) { x = *--X1 - *Y; sm2 += x*x; }
                    *Y = sqrt(sm2*den);
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
