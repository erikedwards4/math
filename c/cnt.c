//Vec2scalar (reduction) operation.
//Gets count of nonzero values for each vector in X along dim.
//This is also the Hamming norm (or L0 norm) of each vector in X.

//For complex case, counts if either real or imaginary part is nonzero.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cnt_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cnt_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cnt_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int cnt_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int cnt_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cnt_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = (float)(*X!=0.0f); }
    }
    else if (L==N)
    {
        size_t cnt = 0u;
        for (size_t l=0u; l<L; ++l, ++X) { cnt += (*X!=0.0f); }
        *Y = cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            size_t cnt;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                cnt = 0u;
                for (size_t l=0u; l<L; ++l, ++X) { cnt += (*X!=0.0f); }
                *Y = cnt;
            }
        }
        else if (G==1u)
        {
            for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y = (float)(*X!=0.0f); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += (float)(*X!=0.0f); }
            }
        }
        else
        {
            size_t cnt;
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    cnt = 0u;
                    for (size_t l=0u; l<L; ++l, X+=K) { cnt += (*X!=0.0f); }
                    *Y = cnt;
                }
            }
        }
    }

    return 0;
}


int cnt_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cnt_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = (double)(*X!=0.0); }
    }
    else if (L==N)
    {
        size_t cnt = 0u;
        for (size_t l=0u; l<L; ++l, ++X) { cnt += (*X!=0.0); }
        *Y = cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            size_t cnt;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                cnt = 0u;
                for (size_t l=0u; l<L; ++l, ++X) { cnt += (*X!=0.0); }
                *Y = cnt;
            }
        }
        else if (G==1u)
        {
            for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y = (double)(*X!=0.0); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, ++X, ++Y) { *Y += (double)(*X!=0.0); }
            }
        }
        else
        {
            size_t cnt;
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    cnt = 0u;
                    for (size_t l=0u; l<L; ++l, X+=K) { cnt += (*X!=0.0); }
                    *Y = cnt;
                }
            }
        }
    }

    return 0;
}


int cnt_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cnt_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, X+=2, ++Y) { *Y = (float)(*X!=0.0f || *(X+1)!=0.0f); }
    }
    else if (L==N)
    {
        size_t cnt = 0u;
        for (size_t l=0u; l<L; ++l, X+=2) { cnt += (*X!=0.0f || *(X+1)!=0.0f); }
        *Y = cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            size_t cnt;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                cnt = 0u;
                for (size_t l=0u; l<L; ++l, X+=2) { cnt += (*X!=0.0f || *(X+1)!=0.0f); }
                *Y = cnt;
            }
        }
        else if (G==1u)
        {
            for (size_t v=0u; v<V; ++v, X+=2, ++Y) { *Y = (float)(*X!=0.0f || *(X+1)!=0.0f); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, X+=2, ++Y) { *Y += (float)(*X!=0.0f || *(X+1)!=0.0f); }
            }
        }
        else
        {
            size_t cnt;
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, ++Y)
                {
                    cnt = 0u;
                    for (size_t l=0u; l<L; ++l, X+=2u*K) { cnt += (*X!=0.0f || *(X+1)!=0.0f); }
                    *Y = cnt;
                }
            }
        }
    }

    return 0;
}


int cnt_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cnt_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, X+=2, ++Y) { *Y = (double)(*X!=0.0 || *(X+1)!=0.0); }
    }
    else if (L==N)
    {
        size_t cnt = 0u;
        for (size_t l=0u; l<L; ++l, X+=2) { cnt += (*X!=0.0 || *(X+1)!=0.0); }
        *Y = cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            size_t cnt;
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                cnt = 0u;
                for (size_t l=0u; l<L; ++l, X+=2) { cnt += (*X!=0.0 || *(X+1)!=0.0); }
                *Y = cnt;
            }
        }
        else if (G==1u)
        {
            for (size_t v=0u; v<V; ++v, X+=2, ++Y) { *Y = (double)(*X!=0.0 || *(X+1)!=0.0); }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, X+=2, ++Y) { *Y += (double)(*X!=0.0 || *(X+1)!=0.0); }
            }
        }
        else
        {
            size_t cnt;
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, ++Y)
                {
                    cnt = 0u;
                    for (size_t l=0u; l<L; ++l, X+=2u*K) { cnt += (*X!=0.0 || *(X+1)!=0.0); }
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
