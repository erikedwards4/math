//Vec2scalar (reduction) operation.
//Gets L0 (Hamming) norm for each vector in X along dim.
//This is just the count of nonzero values for each vector in X.
//This is identical to the cnt function,
//except that I use a thresh of 2*EPS (in case of numerical rounding errors).

//For complex case, counts if either real or imaginary part is nonzero.

#include <stdio.h>
#include <float.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int norm0_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int norm0_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int norm0_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int norm0_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int norm0_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm0_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float pthresh = 2.0f * FLT_EPSILON, nthresh = -pthresh;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = (float)(*X>pthresh || *X<nthresh); }
    }
    else if (L==N)
    {
        size_t cnt = 0;
        for (size_t l=0u; l<L; ++l, ++X) { cnt += (*X>pthresh || *X<nthresh); }
        *Y = cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            size_t cnt;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                cnt = 0;
                for (size_t l=0u; l<L; ++l, ++X) { cnt += (*X>pthresh || *X<nthresh); }
                *Y = cnt;
            }
        }
        else if (G==1)
        {
            for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y = (float)(*X>pthresh || *X<nthresh); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += (float)(*X>pthresh || *X<nthresh); }
            }
        }
        else
        {
            size_t cnt;
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
                {
                    cnt = 0;
                    for (size_t l=0u; l<L; ++l, X+=K) { cnt += (*X>pthresh || *X<nthresh); }
                    *Y = cnt;
                }
            }
        }
    }

    return 0;
}


int norm0_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm0_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double pthresh = 2.0 * DBL_EPSILON, nthresh = -pthresh;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = (double)(*X>pthresh || *X<nthresh); }
    }
    else if (L==N)
    {
        size_t cnt = 0;
        for (size_t l=0u; l<L; ++l, ++X) { cnt += (*X>pthresh || *X<nthresh); }
        *Y = cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            size_t cnt;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                cnt = 0;
                for (size_t l=0u; l<L; ++l, ++X) { cnt += (*X>pthresh || *X<nthresh); }
                *Y = cnt;
            }
        }
        else if (G==1)
        {
            for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y = (double)(*X>pthresh || *X<nthresh); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += (double)(*X>pthresh || *X<nthresh); }
            }
        }
        else
        {
            size_t cnt;
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
                {
                    cnt = 0;
                    for (size_t l=0u; l<L; ++l, X+=K) { cnt += (*X>pthresh || *X<nthresh); }
                    *Y = cnt;
                }
            }
        }
    }

    return 0;
}


int norm0_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm0_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float pthresh = 2.0f * FLT_EPSILON, nthresh = -pthresh;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, X+=2, ++Y) { *Y = (float)(*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
    }
    else if (L==N)
    {
        size_t cnt = 0;
        for (size_t l=0u; l<L; ++l, X+=2) { cnt += (*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
        *Y = cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            size_t cnt;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                cnt = 0;
                for (size_t l=0u; l<L; ++l, X+=2) { cnt += (*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
                *Y = cnt;
            }
        }
        else if (G==1)
        {
            for (size_t v=0u; v<V; ++v, X+=2, ++Y) { *Y = (float)(*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, X+=2, ++Y) { *Y += (float)(*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
            }
        }
        else
        {
            size_t cnt;
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2, ++Y)
                {
                    cnt = 0;
                    for (size_t l=0u; l<L; ++l, X+=2*K) { cnt += (*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
                    *Y = cnt;
                }
            }
        }
    }

    return 0;
}


int norm0_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in norm0_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double pthresh = 2.0 * DBL_EPSILON, nthresh = -pthresh;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, X+=2, ++Y) { *Y = (double)(*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
    }
    else if (L==N)
    {
        size_t cnt = 0;
        for (size_t l=0u; l<L; ++l, X+=2) { cnt += (*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
        *Y = cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            size_t cnt;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                cnt = 0;
                for (size_t l=0u; l<L; ++l, X+=2) { cnt += (*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
                *Y = cnt;
            }
        }
        else if (G==1)
        {
            for (size_t v=0u; v<V; ++v, X+=2, ++Y) { *Y = (double)(*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, X+=2, ++Y) { *Y += (double)(*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
            }
        }
        else
        {
            size_t cnt;
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2, ++Y)
                {
                    cnt = 0;
                    for (size_t l=0u; l<L; ++l, X+=2*K) { cnt += (*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
                    *Y = cnt;
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
