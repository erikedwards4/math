//Vec2scalar (reduction) operation.
//Gets maximum of absolute values for each vector in X along dim.
//This is also the Inf-norm (or max-norm) of each vector.

//For complex case, finds max absolute value and outputs the complex number.
//For complex case, this is |Xr|+|Xi|; see max for the usual sqrt(Xr*Xr+Xl*Xi).
//This according to BLAS standard, but is not actually the Inf-norm.


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int amax_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in amax_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = fabsf(*X); }
    }
    else if (L==N)
    {
        size_t l = cblas_isamax((int)L,X,1);
        *Y = fabsf(*(X+l));
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++Y)
            {
                size_t l = cblas_isamax((int)L,X,1);
                X += l; *Y = fabsf(*X); X += L-l;
            }
        }
        else if (G==1u)
        {
            float *mxs;
            if (!(mxs=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in amax_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=L; l>0u; --l, Y-=V, mxs-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y, ++mxs)
                {
                    if (*X**X>*mxs) { *Y = fabsf(*X); *mxs = *X**X; }
                }
            }
            free(mxs);
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++Y)
                {
                    size_t l = cblas_isamax((int)L,X,(int)K);
                    X += l*K; *Y = fabsf(*X);
                    X -= (int)((L-l)*K)-1;
                }
            }
        }
    }

    return 0;
}


int amax_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in amax_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = fabs(*X); }
    }
    else if (L==N)
    {
        size_t l = cblas_idamax((int)L,X,1);
        *Y = fabs(*(X+l));
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++Y)
            {
                size_t l = cblas_idamax((int)L,X,1);
                X += l; *Y = fabs(*X); X += L-l;
            }
        }
        else if (G==1u)
        {
            double *mxs;
            if (!(mxs=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in amax_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=L; l>0u; --l, Y-=V, mxs-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y, ++mxs)
                {
                    if (*X**X>*mxs) { *Y = fabs(*X); *mxs = *X**X; }
                }
            }
            free(mxs);
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++Y)
                {
                    size_t l = cblas_idamax((int)L,X,(int)K);
                    X += l*K; *Y = fabs(*X);
                    X -= (int)((L-l)*K)-1;
                }
            }
        }
    }

    return 0;
}


int amax_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in amax_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (L==1u)
    {
        for (size_t n=N; n>0u; --n, X+=2, ++Y) { *Y = fabsf(*X) + fabsf(*(X+1)); }
    }
    else if (L==N)
    {
        X += 2 * cblas_icamax((int)L,X,1);
        *Y = fabsf(*X) + fabsf(*(X+1));
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
                size_t l = cblas_icamax((int)L,X,1);
                X += 2u*l;
                *Y++ = fabsf(*X) + fabsf(*(X+1));
                X += 2u*(L-l);
            }
        }
        else if (G==1u)
        {
            for (size_t l=L; l>0u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, X+=2, ++Y)
                {
                    float ax = fabsf(*X) + fabsf(*(X+1));
                    if (l==0u || ax>*Y) { *Y = ax;}
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b)
                {
                    size_t l = cblas_icamax((int)L,X,(int)K);
                    X += 2u*l*K;
                    *Y++ = fabsf(*X) + fabsf(*(X+1));
                    X -= 2*((int)((L-l)*K)-1);
                }
            }
        }
    }

    return 0;
}


int amax_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in amax_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (L==1u)
    {
        for (size_t n=N; n>0u; --n, X+=2, ++Y) { *Y = fabs(*X) + fabs(*(X+1)); }
    }
    else if (L==N)
    {
        X += 2 * cblas_izamax((int)L,X,1);
        *Y = fabs(*X) + fabs(*(X+1));
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
                size_t l = cblas_izamax((int)L,X,1);
                X += 2u*l;
                *Y++ = fabs(*X) + fabs(*(X+1));
                X += 2u*(L-l);
            }
        }
        else if (G==1u)
        {
            for (size_t l=L; l>0u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, X+=2, ++Y)
                {
                    double ax = fabs(*X) + fabs(*(X+1));
                    if (l==0u || ax>*Y) { *Y = ax;}
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b)
                {
                    size_t l = cblas_izamax((int)L,X,(int)K);
                    X += 2u*l*K;
                    *Y++ = fabs(*X) + fabs(*(X+1));
                    X -= 2*((int)((L-l)*K)-1);
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
