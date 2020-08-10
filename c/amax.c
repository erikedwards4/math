//Vec2scalar (reduction) operation.
//Gets maximum of absolute values for each vector in X along dim.
//This is also the Inf-norm (or max-norm) of each vector.

//For complex case, finds max absolute value and outputs the complex number.
//For complex case, this is |Xr|+|Xi|; see max for the usual sqrt(Xr*Xr+Xl*Xi).
//This according to BLAS standard, but is not actually the Inf-norm.


#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int amax_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int amax_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int amax_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int amax_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int amax_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in amax_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        size_t l = cblas_isamax((int)L,X,1);
        *Y = *(X+l);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, ++Y)
            {
                size_t l = cblas_isamax((int)L,X,1);
                X += l; *Y = *X; X += L-l;
            }
        }
        else if (G==1)
        {
            float *mxs;
            if (!(mxs=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in amax_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0; l<L; ++l, Y-=V, mxs-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y, ++mxs)
                {
                    if (*X**X>*mxs) { *Y = *X; *mxs = *X**X; }
                }
            }
            free(mxs);
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++Y)
                {
                    size_t l = cblas_isamax((int)L,X,(int)K);
                    X += l*K; *Y = *X;
                    X -= (int)((L-l)*K)-1;
                }
            }
        }
    }

    return 0;
}


int amax_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in amax_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        size_t l = cblas_idamax((int)L,X,1);
        *Y = *(X+l);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, ++Y)
            {
                size_t l = cblas_idamax((int)L,X,1);
                X += l; *Y = *X; X += L-l;
            }
        }
        else if (G==1)
        {
            double *mxs;
            if (!(mxs=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in amax_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0; l<L; ++l, Y-=V, mxs-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y, ++mxs)
                {
                    if (*X**X>*mxs) { *Y = *X; *mxs = *X**X; }
                }
            }
            free(mxs);
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++Y)
                {
                    size_t l = cblas_idamax((int)L,X,(int)K);
                    X += l*K; *Y = *X;
                    X -= (int)((L-l)*K)-1;
                }
            }
        }
    }

    return 0;
}


int amax_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in amax_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (L==1)
    {
        for (size_t n=0; n<N; ++n, X+=2, ++Y) { *Y = fabsf(*X) + fabsf(*(X+1)); }
    }
    else if (L==N)
    {
        X += 2 * cblas_icamax((int)L,X,1);
        *Y = fabsf(*X) + fabsf(*(X+1));
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                size_t l = cblas_icamax((int)L,X,1);
                X += 2*l;
                *Y++ = fabsf(*X) + fabsf(*(X+1));
                X += 2*(L-l);
            }
        }
        else if (G==1)
        {
            for (size_t l=0; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, X+=2, ++Y)
                {
                    float ax = fabsf(*X) + fabsf(*(X+1));
                    if (l==0 || ax>*Y) { *Y = ax;}
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b)
                {
                    size_t l = cblas_icamax((int)L,X,(int)K);
                    X += 2*l*K;
                    *Y++ = fabsf(*X) + fabsf(*(X+1));
                    X -= 2*((int)((L-l)*K)-1);
                }
            }
        }
    }

    return 0;
}


int amax_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in amax_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (L==1)
    {
        for (size_t n=0; n<N; ++n, X+=2, ++Y) { *Y = fabs(*X) + fabs(*(X+1)); }
    }
    else if (L==N)
    {
        X += 2 * cblas_izamax((int)L,X,1);
        *Y = fabs(*X) + fabs(*(X+1));
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                size_t l = cblas_izamax((int)L,X,1);
                X += 2*l;
                *Y++ = fabs(*X) + fabs(*(X+1));
                X += 2*(L-l);
            }
        }
        else if (G==1)
        {
            for (size_t l=0; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, X+=2, ++Y)
                {
                    double ax = fabs(*X) + fabs(*(X+1));
                    if (l==0 || ax>*Y) { *Y = ax;}
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b)
                {
                    size_t l = cblas_izamax((int)L,X,(int)K);
                    X += 2*l*K;
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
