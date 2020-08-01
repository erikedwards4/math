//Vec2scalar (reduction) operation.
//Gets geometric mean for each vector in X along dim.
//This is the Lth root of the prod for each vector.
//This is also the exp of the mean of logs for each vector.

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int geomean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int geomean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int geomean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int geomean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int geomean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in geomean_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y = 0.0f;
        for (size_t l=0; l<L; ++l, ++X) { *Y += logf(*X); }
        *Y = expf(*Y/L);
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
                *Y = 0.0f;
                for (size_t l=0; l<L; ++l, ++X) { *Y += logf(*X); }
                *Y = expf(*Y/L);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = logf(*X); }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += logf(*X); }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = expf(*Y/L); }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    *Y = 0.0f;
                    for (size_t l=0; l<L; ++l, X+=K) { *Y += logf(*X); }
                    *Y = expf(*Y/L);
                }
            }
        }
    }

    return 0;
}


int geomean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in geomean_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y = 0.0;
        for (size_t l=0; l<L; ++l, ++X) { *Y += log(*X); }
        *Y = exp(*Y/L);
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
                *Y = 0.0;
                for (size_t l=0; l<L; ++l, ++X) { *Y += log(*X); }
                *Y = exp(*Y/L);
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = log(*X); }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += log(*X); }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = exp(*Y/L); }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    *Y = 0.0;
                    for (size_t l=0; l<L; ++l, X+=K) { *Y += log(*X); }
                    *Y = exp(*Y/L);
                }
            }
        }
    }

    return 0;
}


int geomean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in geomean_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    _Complex float y;

    if (N==0) {}
    else if (L==1) { cblas_ccopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y++ = 0.0f; *Y-- = 0.0f;
        for (size_t l=0; l<L; ++l, X+=2)
        {
            y = clogf(*X + 1.0if**(X+1));
            *Y++ += *(float *)&y; *Y-- += *((float *)&y+1);
        }
        y = cexpf(*Y/L + 1.0if**(Y+1)/L);
        *Y++ = *(float *)&y; *Y = *((float *)&y+1);
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
                *Y++ = 0.0f; *Y-- = 0.0f;
                for (size_t l=0; l<L; ++l, X+=2)
                {
                    y = clogf(*X + 1.0if**(X+1));
                    *Y++ += *(float *)&y; *Y-- += *((float *)&y+1);
                }
                y = cexpf(*Y/L + 1.0if**(Y+1)/L);
                *Y++ = *(float *)&y; *Y++ = *((float *)&y+1);
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2)
                {
                    *Y++ = 0.0f; *Y-- = 0.0f;
                    for (size_t l=0; l<L; ++l, X+=2*K)
                    {
                        y = clogf(*X + 1.0if**(X+1));
                        *Y++ += *(float *)&y; *Y-- += *((float *)&y+1);
                    }
                    y = cexpf(*Y/L + 1.0if**(Y+1)/L);
                    *Y++ = *(float *)&y; *Y++ = *((float *)&y+1);
                }
            }
        }
    }

    return 0;
}


int geomean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in geomean_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    _Complex double y;

    if (N==0) {}
    else if (L==1) { cblas_zcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y++ = 0.0; *Y-- = 0.0;
        for (size_t l=0; l<L; ++l, X+=2)
        {
            y = clog(*X + 1.0i**(X+1));
            *Y++ += *(double *)&y; *Y-- += *((double *)&y+1);
        }
        y = cexp(*Y/L + 1.0i**(Y+1)/L);
        *Y++ = *(double *)&y; *Y = *((double *)&y+1);
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
                *Y++ = 0.0; *Y-- = 0.0;
                for (size_t l=0; l<L; ++l, X+=2)
                {
                    y = clog(*X + 1.0i**(X+1));
                    *Y++ += *(double *)&y; *Y-- += *((double *)&y+1);
                }
                y = cexp(*Y/L + 1.0i**(Y+1)/L);
                *Y++ = *(double *)&y; *Y++ = *((double *)&y+1);
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2)
                {
                    *Y++ = 0.0; *Y-- = 0.0;
                    for (size_t l=0; l<L; ++l, X+=2*K)
                    {
                        y = clog(*X + 1.0i**(X+1));
                        *Y++ += *(double *)&y; *Y-- += *((double *)&y+1);
                    }
                    y = cexp(*Y/L + 1.0i**(Y+1)/L);
                    *Y++ = *(double *)&y; *Y++ = *((double *)&y+1);
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
