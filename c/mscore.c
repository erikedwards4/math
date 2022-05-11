//M-scores each vector in X along dim.
//Subtracts median and divides by the MAD (median absolute deviation from the median).

#include <stdio.h>
//#include <lapacke.h>
#include "codee_math.h"
#include "extremum.c"
#include "kselect.c"


#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int mscore_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in mscore_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in mscore_s: L (vec length) must be > 1\n"); return 1; }
    const size_t i2 = L/2u;
    float med, mad;
    
    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        Y -= L; X -= L;
        if (LAPACKE_slasrt_work('I',(int)L,Y)) { fprintf(stderr,"error in mscore_s: problem with LAPACKE function\n"); }
        Y += i2;
        med = (L%2u) ? *Y : 0.5f*(*Y + *(Y-1));
        Y -= i2;
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = fabsf(*X-med); }
        Y -= L;
        if (LAPACKE_slasrt_work('I',(int)L,Y)) { fprintf(stderr,"error in mscore_s: problem with LAPACKE function\n"); }
        Y += i2;
        mad = (L%2u) ? *Y : 0.5f*(*Y + *(Y-1));
        Y -= i2;
        for (size_t l=L; l>0u; --l, ++Y) { *Y /= mad; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;
        float *X1;
        if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in mscore_s: problem with malloc. "); perror("malloc"); return 1; }

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_s: problem with LAPACKE function\n"); }
                X1 += i2;
                med = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=L; l>0u; --l, ++X, ++X1, ++Y) { *Y = *X-med; *X1 = fabsf(*Y); }
                X1 -= L; Y -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_s: problem with LAPACKE function\n"); }
                X1 += i2;
                mad = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=L; l>0u; --l, ++Y) { *Y /= mad; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_s: problem with LAPACKE function\n"); }
                    X1 += i2;
                    med = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=L; l>0u; --l, X+=K, Y+=K, ++X1) { *Y = *X-med; *X1 = fabsf(*Y); }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_s: problem with LAPACKE function\n"); }
                    X1 += i2;
                    mad = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=L; l>0u; --l, Y+=K) { *Y /= mad; }
                }
            }
        }
        free(X1);
    }

    return 0;
}


int mscore_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in mscore_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in mscore_d: L (vec length) must be > 1\n"); return 1; }
    const size_t i2 = L/2u;
    double med, mad;

    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = *X; }
        Y -= L; X -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,Y)) { fprintf(stderr,"error in mscore_d: problem with LAPACKE function\n"); }
        Y += i2;
        med = (L%2u) ? *Y : 0.5*(*Y + *(Y-1));
        Y -= i2;
        for (size_t l=L; l>0u; --l, ++X, ++Y) { *Y = fabs(*X-med); }
        Y -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,Y)) { fprintf(stderr,"error in mscore_d: problem with LAPACKE function\n"); }
        Y += i2;
        mad = (L%2u) ? *Y : 0.5*(*Y + *(Y-1));
        Y -= i2;
        for (size_t l=L; l>0u; --l, ++Y) { *Y /= mad; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;
        double *X1;
        if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in mscore_d: problem with malloc. "); perror("malloc"); return 1; }

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_d: problem with LAPACKE function\n"); }
                X1 += i2;
                med = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=L; l>0u; --l, ++X, ++X1, ++Y) { *Y = *X-med; *X1 = fabs(*Y); }
                X1 -= L; Y -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_d: problem with LAPACKE function\n"); }
                X1 += i2;
                mad = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=L; l>0u; --l, ++Y) { *Y /= mad; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y-=K*L-1u)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_d: problem with LAPACKE function\n"); }
                    X1 += i2;
                    med = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=L; l>0u; --l, X+=K, Y+=K, ++X1) { *Y = *X-med; *X1 = fabs(*Y); }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_d: problem with LAPACKE function\n"); }
                    X1 += i2;
                    mad = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=L; l>0u; --l, Y+=K) { *Y /= mad; }
                }
            }
        }
        free(X1);
    }

    return 0;
}


int mscore_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in mscore_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in mscore_inplace_s: L (vec length) must be > 1\n"); return 1; }
    const size_t i2 = L/2u;
    float med, mad;

    float *X1;
    if (!(X1=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in mscore_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L; X -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_inplace_s: problem with LAPACKE function\n"); }
        X1 += i2;
        med = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
        X1 -= i2;
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X -= med; *X1 = fabsf(*X); }
        X1 -= L; X -= L;
        if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_inplace_s: problem with LAPACKE function\n"); }
        X1 += i2;
        mad = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
        X1 -= i2;
        for (size_t l=L; l>0u; --l, ++X) { *X /= mad; }
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
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_inplace_s: problem with LAPACKE function\n"); }
                X1 += i2;
                med = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X -= med; *X1 = fabsf(*X); }
                X1 -= L; X -= L;
                if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_inplace_s: problem with LAPACKE function\n"); }
                X1 += i2;
                mad = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=L; l>0u; --l, ++X) { *X /= mad; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L; X -= K*L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_inplace_s: problem with LAPACKE function\n"); }
                    X1 += i2;
                    med = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X -= med; *X1 = fabsf(*X); }
                    X1 -= L;
                    if (LAPACKE_slasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_inplace_s: problem with LAPACKE function\n"); }
                    X1 += i2;
                    mad = (L%2u) ? *X1 : 0.5f*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=L; l>0u; --l) { X-=K; *X /= mad; }
                }
            }
        }
    }

    free(X1);
    return 0;
}


int mscore_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in mscore_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L<2u) { fprintf(stderr,"error in mscore_inplace_d: L (vec length) must be > 1\n"); return 1; }
    const size_t i2 = L/2u;
    double med, mad;

    double *X1;
    if (!(X1=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in mscore_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N==0u) {}
    else if (L==N)
    {
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
        X1 -= L; X -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_inplace_d: problem with LAPACKE function\n"); }
        X1 += i2;
        med = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
        X1 -= i2;
        for (size_t l=L; l>0u; --l, ++X, ++X1) { *X -= med; *X1 = fabs(*X); }
        X1 -= L; X -= L;
        if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_inplace_d: problem with LAPACKE function\n"); }
        X1 += i2;
        mad = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
        X1 -= i2;
        for (size_t l=L; l>0u; --l, ++X) { *X /= mad; }
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
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= L; X -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_inplace_d: problem with LAPACKE function\n"); }
                X1 += i2;
                med = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=L; l>0u; --l, ++X, ++X1) { *X -= med; *X1 = fabs(*X); }
                X1 -= L; X -= L;
                if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_inplace_d: problem with LAPACKE function\n"); }
                X1 += i2;
                mad = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
                X1 -= i2;
                for (size_t l=L; l>0u; --l, ++X) { *X /= mad; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X)
                {
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L; X -= K*L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_inplace_d: problem with LAPACKE function\n"); }
                    X1 += i2;
                    med = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X -= med; *X1 = fabs(*X); }
                    X1 -= L;
                    if (LAPACKE_dlasrt_work('I',(int)L,X1)) { fprintf(stderr,"error in mscore_inplace_d: problem with LAPACKE function\n"); }
                    X1 += i2;
                    mad = (L%2u) ? *X1 : 0.5*(*X1 + *(X1-1));
                    X1 -= i2;
                    for (size_t l=L; l>0u; --l) { X-=K; *X /= mad; }
                }
            }
        }
    }

    free(X1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
