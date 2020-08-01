//Vec2scalar (reduction) operation.
//Gets sum of absolute-values (L1-norm) for each vector in X along dim.
//For complex case, output is real and is sum(|Xr|+|Xi|), not sum(|X|). See norm1 for the later.

//Thus, this code is identical to norm1 for real-valued case, but not for complex.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int asum_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int asum_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int asum_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int asum_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int asum_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in asum_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = fabsf(*X); }
    }
    else if (L==N)
    {
        *Y = cblas_sasum((int)L,X,1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L, ++Y)
            {
               *Y = cblas_sasum((int)L,X,1);
            }
        }
        else if (G==1)
        {
            // for (size_t v=0; v<V; ++v, ++X, ++Y)
            // {
            //    *Y = cblas_sasum((int)L,X,(int)V);
            // }
            const float z = 0.0f;
            cblas_scopy((int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += fabsf(*X); }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                   *Y = cblas_sasum((int)L,X,(int)K);
                }
            }
        }
    }

    return 0;
}


int asum_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in asum_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = fabs(*X); }
    }
    else if (L==N)
    {
        *Y = cblas_dasum((int)L,X,1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=L, ++Y)
            {
               *Y = cblas_dasum((int)L,X,1);
            }
        }
        else if (G==1)
        {
            const double z = 0.0;
            cblas_dcopy((int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += fabs(*X); }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, ++X, ++Y)
                {
                   *Y = cblas_dasum((int)L,X,(int)K);
                }
            }
        }
    }

    return 0;
}


int asum_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in asum_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, X+=2, ++Y) { *Y = fabsf(*X) + fabsf(*(X+1)); }
    }
    else if (L==N)
    {
        *Y = cblas_scasum((int)L,X,1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=2*L, ++Y)
            {
               *Y = cblas_scasum((int)L,X,1);
            }
        }
        else if (G==1)
        {
            const float z = 0.0f;
            cblas_scopy((int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, X+=2, ++Y) { *Y += fabsf(*X) + fabsf(*(X+1)); }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X+=2, ++Y)
                {
                   *Y = cblas_scasum((int)L,X,(int)K);
                }
            }
        }
    }

    return 0;
}


int asum_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in asum_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, X+=2, ++Y) { *Y = fabs(*X) + fabs(*(X+1)); }
    }
    else if (L==N)
    {
        *Y = cblas_dzasum((int)L,X,1);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, X+=2*L, ++Y)
            {
               *Y = cblas_dzasum((int)L,X,1);
            }
        }
        else if (G==1)
        {
            const double z = 0.0;
            cblas_dcopy((int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, X+=2, ++Y) { *Y += fabs(*X) + fabs(*(X+1)); }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X+=2, ++Y)
                {
                   *Y = cblas_dzasum((int)L,X,(int)K);
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
