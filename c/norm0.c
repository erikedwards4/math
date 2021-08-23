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

int norm0_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int norm0_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int norm0_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int norm0_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int norm0_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in norm0_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float pthresh = 2.0f * FLT_EPSILON, nthresh = -pthresh;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = (float)(*X>pthresh || *X<nthresh); }
    }
    else if (L==N)
    {
        int cnt = 0;
        for (size_t l=L; l>0u; --l, ++X) { cnt += (*X>pthresh || *X<nthresh); }
        *Y = (float)cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            int cnt;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                cnt = 0;
                for (size_t l=L; l>0u; --l, ++X) { cnt += (*X>pthresh || *X<nthresh); }
                *Y = (float)cnt;
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = (float)(*X>pthresh || *X<nthresh); }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += (float)(*X>pthresh || *X<nthresh); }
            }
        }
        else
        {
            int cnt;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    cnt = 0;
                    for (size_t l=L; l>0u; --l, X+=K) { cnt += (*X>pthresh || *X<nthresh); }
                    *Y = (float)cnt;
                }
            }
        }
    }

    return 0;
}


int norm0_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in norm0_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double pthresh = 2.0 * DBL_EPSILON, nthresh = -pthresh;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = (double)(*X>pthresh || *X<nthresh); }
    }
    else if (L==N)
    {
        int cnt = 0;
        for (size_t l=L; l>0u; --l, ++X) { cnt += (*X>pthresh || *X<nthresh); }
        *Y = (double)cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            int cnt;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                cnt = 0;
                for (size_t l=L; l>0u; --l, ++X) { cnt += (*X>pthresh || *X<nthresh); }
                *Y = (double)cnt;
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = (double)(*X>pthresh || *X<nthresh); }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += (double)(*X>pthresh || *X<nthresh); }
            }
        }
        else
        {
            int cnt;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    cnt = 0;
                    for (size_t l=L; l>0u; --l, X+=K) { cnt += (*X>pthresh || *X<nthresh); }
                    *Y = (double)cnt;
                }
            }
        }
    }

    return 0;
}


int norm0_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in norm0_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float pthresh = 2.0f * FLT_EPSILON, nthresh = -pthresh;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, X+=2, ++Y) { *Y = (float)(*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
    }
    else if (L==N)
    {
        int cnt = 0;
        for (size_t l=L; l>0u; --l, X+=2) { cnt += (*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
        *Y = (float)cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            int cnt;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                cnt = 0;
                for (size_t l=L; l>0u; --l, X+=2) { cnt += (*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
                *Y = (float)cnt;
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y = (float)(*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y += (float)(*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
            }
        }
        else
        {
            int cnt;
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    cnt = 0;
                    for (size_t l=L; l>0u; --l, X+=2u*K) { cnt += (*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
                    *Y = (float)cnt;
                }
            }
        }
    }

    return 0;
}


int norm0_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in norm0_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double pthresh = 2.0 * DBL_EPSILON, nthresh = -pthresh;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, X+=2, ++Y) { *Y = (double)(*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
    }
    else if (L==N)
    {
        int cnt = 0;
        for (size_t l=L; l>0u; --l, X+=2) { cnt += (*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
        *Y = (double)cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            int cnt;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                cnt = 0;
                for (size_t l=L; l>0u; --l, X+=2) { cnt += (*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
                *Y = (double)cnt;
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y = (double)(*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y += (double)(*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
            }
        }
        else
        {
            int cnt;
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    cnt = 0;
                    for (size_t l=L; l>0u; --l, X+=2u*K) { cnt += (*X>pthresh || *X<nthresh || *(X+1)>pthresh || *(X+1)<nthresh); }
                    *Y = (double)cnt;
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
