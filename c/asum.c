//Vec2scalar (reduction) operation.
//Gets sum of absolute-values (L1-norm) for each vector in X along dim.
//For complex case, output is real and is sum(|Xr|+|Xi|), not sum(|X|). See norm1 for the later.

//Thus, this code is identical to norm1 for real-valued case, but not for complex.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int asum_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int asum_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int asum_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int asum_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int asum_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in asum_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = fabsf(*X); }
    }
    else if (L==N)
    {
        if (L<20000u)
        {
            float sm = 0.0f;
            for (size_t l=0u; l<L; ++l, ++X) { sm += fabsf(*X); }
            *Y = sm;
        }
        else { *Y = cblas_sasum((int)L,X,1); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (L<25u)
            {
                float sm;
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<L; ++l, ++X) { sm += fabsf(*X); }
                    *Y = sm;
                }
            }
            else
            {
                for (size_t v=V; v>0u; --v, X+=L, ++Y)
                {
                    *Y = cblas_sasum((int)L,X,1);
                }
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = fabsf(*X); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += fabsf(*X); }
            }
        }
        else
        {
            float sm;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<L; ++l, X+=K) { sm += fabsf(*X); }
                    *Y = sm;
                }
            }
        }
    }

    return 0;
}


int asum_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in asum_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = fabs(*X); }
    }
    else if (L==N)
    {
        if (L<20000u)
        {
            double sm = 0.0;
            for (size_t l=0u; l<L; ++l, ++X) { sm += fabs(*X); }
            *Y = sm;
        }
        else { *Y = cblas_dasum((int)L,X,1); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (L<25u)
            {
                double sm;
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<L; ++l, ++X) { sm += fabs(*X); }
                    *Y = sm;
                }
            }
            else
            {
                for (size_t v=V; v>0u; --v, X+=L, ++Y)
                {
                    *Y = cblas_dasum((int)L,X,1);
                }
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = fabs(*X); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += fabs(*X); }
            }
        }
        else
        {
            double sm;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<L; ++l, X+=K) { sm += fabs(*X); }
                    *Y = sm;
                }
            }
        }
    }

    return 0;
}


int asum_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in asum_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, X+=2, ++Y) { *Y = fabsf(*X) + fabsf(*(X+1)); }
    }
    else if (L==N)
    {
        if (L<5000u)
        {
            float sm = 0.0f;
            for (size_t l=0u; l<L; ++l, X+=2) { sm += fabsf(*X) + fabsf(*(X+1)); }
            *Y = sm;
        }
        else { *Y = cblas_scasum((int)L,X,1); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=2u*L, ++Y)
            {
               *Y = cblas_scasum((int)L,X,1);
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y = fabsf(*X) + fabsf(*(X+1)); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y += fabsf(*X) + fabsf(*(X+1)); }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X+=2, ++Y)
                {
                   *Y = cblas_scasum((int)L,X,(int)K);
                }
            }
        }
    }

    return 0;
}


int asum_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in asum_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, X+=2, ++Y) { *Y = fabs(*X) + fabs(*(X+1)); }
    }
    else if (L==N)
    {
        if (L<5000u)
        {
            double sm = 0.0;
            for (size_t l=0u; l<L; ++l, X+=2) { sm += fabs(*X) + fabs(*(X+1)); }
            *Y = sm;
        }
        else { *Y = cblas_dzasum((int)L,X,1); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=2u*L, ++Y)
            {
               *Y = cblas_dzasum((int)L,X,1);
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y = fabs(*X) + fabs(*(X+1)); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y += fabs(*X) + fabs(*(X+1)); }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X+=2, ++Y)
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
