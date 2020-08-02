//Vec2scalar (reduction) operation.
//Gets geometric standard deviation for each vector in X along dim.
//This is the exp of the std of logs for each vector.

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int geostd_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int geostd_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int geostd_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int geostd_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int geostd_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in geostd_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f/L, den2 = den;
    float x, mn, sm2;

    if (N==0) {}
    else if (L==1)
    {
        const float o = 1.0f;
        cblas_scopy((int)N,&o,0,Y,1);
    }
    else if (L==N)
    {
        mn = sm2 = 0.0f;
        if (L<150)
        {
            for (size_t l=0; l<L; ++l, ++X) { mn += logf(*X); }
            mn *= den;
            for (size_t l=0; l<L; ++l) { x = logf(*--X) - mn; sm2 += x*x; }
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in geostd_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t l=0; l<L; ++l, ++X1, ++X) { *X1 = logf(*X); mn += *X1; }
            mn *= den;
            for (size_t l=0; l<L; ++l) { x = *--X1 - mn; sm2 += x*x; }
            free(X1);
        }
        *Y = expf(sm2*den2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in geostd_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t v=0; v<V; ++v, ++Y)
            {
                mn = sm2 = 0.0f;
                for (size_t l=0; l<L; ++l, ++X1, ++X) { *X1 = logf(*X); mn += *X1; }
                mn *= den;
                for (size_t l=0; l<L; ++l) { x = *--X1 - mn; sm2 += x*x; }
                *Y = expf(sm2*den2);
            }
            free(X1);
        }
        else
        {
            float *X1;
            if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in geostd_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    mn = sm2 = 0.0f;
                    for (size_t l=0; l<L; ++l, ++X1, X+=K) { *X1 = logf(*X); mn += *X1; }
                    mn *= den;
                    for (size_t l=0; l<L; ++l) { x = *--X1 - mn; sm2 += x*x; }
                    *Y = expf(sm2*den2);
                }
            }
            free(X1);
        }
    }

    return 0;
}


int geostd_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in geostd_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0/L, den2 = den;
    double x, mn, sm2;

    if (N==0) {}
    else if (L==1)
    {
        const double o = 1.0;
        cblas_dcopy((int)N,&o,0,Y,1);
    }
    else if (L==N)
    {
        mn = sm2 = 0.0;
        if (L<150)
        {
            for (size_t l=0; l<L; ++l, ++X) { mn += log(*X); }
            mn *= den;
            for (size_t l=0; l<L; ++l) { x = log(*--X) - mn; sm2 += x*x; }
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in geostd_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t l=0; l<L; ++l, ++X1, ++X) { *X1 = log(*X); mn += *X1; }
            mn *= den;
            for (size_t l=0; l<L; ++l) { x = *--X1 - mn; sm2 += x*x; }
            free(X1);
        }
        *Y = exp(sm2*den2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in geostd_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t v=0; v<V; ++v, ++Y)
            {
                mn = sm2 = 0.0;
                for (size_t l=0; l<L; ++l, ++X1, ++X) { *X1 = log(*X); mn += *X1; }
                mn *= den;
                for (size_t l=0; l<L; ++l) { x = *--X1 - mn; sm2 += x*x; }
                *Y = exp(sm2*den2);
            }
            free(X1);
        }
        else
        {
            double *X1;
            if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in geostd_d: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    mn = sm2 = 0.0;
                    for (size_t l=0; l<L; ++l, ++X1, X+=K) { *X1 = log(*X); mn += *X1; }
                    mn *= den;
                    for (size_t l=0; l<L; ++l) { x = *--X1 - mn; sm2 += x*x; }
                    *Y = exp(sm2*den2);
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
