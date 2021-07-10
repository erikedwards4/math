//M-scores each vector in X along dim.
//Subtracts median and divides by the MAD (median absolute deviation from the median).
//This operates in-place.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>


#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mscore_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mscore_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int mscore_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mscore_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (L<2) { fprintf(stderr,"error in mscore_s: L (vec length) must be > 1\n"); return 1; }
    const size_t i2 = L/2;
    float med, mad;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in mscore_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L; X -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_s: problem with LAPACKE function\n"); }
        X1 += i2;
        med = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
        X1 -= i2;
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X -= med; *X1 = fabsf(*X); }
        X1 -= L; X -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_s: problem with LAPACKE function\n"); }
        X1 += i2;
        mad = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
        X1 -= i2;
        for (size_t l=0u; l<L; ++l, ++X) { *X /= mad; }
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
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_s: problem with LAPACKE function\n"); }
                X1 += i2;
                med = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X -= med; *X1 = fabsf(*X); }
                X1 -= L; X -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_s: problem with LAPACKE function\n"); }
                X1 += i2;
                mad = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=0u; l<L; ++l, ++X) { *X /= mad; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, ++X)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L; X -= K*L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_s: problem with LAPACKE function\n"); }
                    X1 += i2;
                    med = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X -= med; *X1 = fabsf(*X); }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_s: problem with LAPACKE function\n"); }
                    X1 += i2;
                    mad = (L%2) ? *X1 : 0.5f*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=0u; l<L; ++l) { X-=K; *X /= mad; }
                }
            }
        }
    }

    return 0;
}


int mscore_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mscore_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (L<2) { fprintf(stderr,"error in mscore_d: L (vec length) must be > 1\n"); return 1; }
    const size_t i2 = L/2;
    double med, mad;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in mscore_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= L; X -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_d: problem with LAPACKE function\n"); }
        X1 += i2;
        med = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
        X1 -= i2;
        for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X -= med; *X1 = fabs(*X); }
        X1 -= L; X -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_d: problem with LAPACKE function\n"); }
        X1 += i2;
        mad = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
        X1 -= i2;
        for (size_t l=0u; l<L; ++l, ++X) { *X /= mad; }
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
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_d: problem with LAPACKE function\n"); }
                X1 += i2;
                med = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X -= med; *X1 = fabs(*X); }
                X1 -= L; X -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_d: problem with LAPACKE function\n"); }
                X1 += i2;
                mad = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=0u; l<L; ++l, ++X) { *X /= mad; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, ++X)
                {
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L; X -= K*L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_d: problem with LAPACKE function\n"); }
                    X1 += i2;
                    med = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X -= med; *X1 = fabs(*X); }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_d: problem with LAPACKE function\n"); }
                    X1 += i2;
                    mad = (L%2) ? *X1 : 0.5*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=0u; l<L; ++l) { X-=K; *X /= mad; }
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
