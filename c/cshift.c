//Circular shift of X by P elements along dim.
//This has no in-place version (no faster version that I can find).

#include <stdio.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cshift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P);
int cshift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P);
int cshift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P);
int cshift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P);


int cshift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in cshift_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0u) {}
    else if (D==0)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        X += L2;
        for (size_t l=0u; l<L1; ++l, ++X, ++Y) { *Y = *X; }
        X -= L;
        for (size_t l=0u; l<L2; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += L2;
            for (size_t v=0u; v<V; ++v, X+=L)
            {
                for (size_t l=0u; l<L1; ++l, ++X, ++Y) { *Y = *X; }
                X -= L;
                for (size_t l=0u; l<L2; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            X += V * L2;
            for (size_t n=0u; n<V*L1; ++n, ++X, ++Y) { *Y = *X; }
            X -= N;
            for (size_t n=0u; n<V*L2; ++n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                X += K * L2;
                for (size_t b=0u; b<B; ++b, X-=K*L2-1u, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<L1; ++l, X+=K, Y+=K) { *Y = *X; }
                    X -= K*L;
                    for (size_t l=0u; l<L2; ++l, X+=K, Y+=K) { *Y = *X; }
                }
            }
        }
    }
    
    return 0;
}


int cshift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in cshift_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0u) {}
    else if (D==0)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        X += L2;
        for (size_t l=0u; l<L1; ++l, ++X, ++Y) { *Y = *X; }
        X -= L;
        for (size_t l=0u; l<L2; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += L2;
            for (size_t v=0u; v<V; ++v, X+=L)
            {
                for (size_t l=0u; l<L1; ++l, ++X, ++Y) { *Y = *X; }
                X -= L;
                for (size_t l=0u; l<L2; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            X += V * L2;
            for (size_t n=0u; n<V*L1; ++n, ++X, ++Y) { *Y = *X; }
            X -= N;
            for (size_t n=0u; n<V*L2; ++n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                X += K * L2;
                for (size_t b=0u; b<B; ++b, X-=K*L2-1u, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<L1; ++l, X+=K, Y+=K) { *Y = *X; }
                    X -= K*L;
                    for (size_t l=0u; l<L2; ++l, X+=K, Y+=K) { *Y = *X; }
                }
            }
        }
    }
    
    return 0;
}


int cshift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in cshift_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0u) {}
    else if (D==0)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        X += 2u*L2;
        for (size_t l=0u; l<2u*L1; ++l, ++X, ++Y) { *Y = *X; }
        X -= 2u*L;
        for (size_t l=0u; l<2u*L2; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += 2u*L2;
            for (size_t v=0u; v<V; ++v, X+=2u*L)
            {
                for (size_t l=0u; l<2u*L1; ++l, ++X, ++Y) { *Y = *X; }
                X -= 2u*L;
                for (size_t l=0u; l<2u*L2; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            X += 2u*V*L2;
            for (size_t n=0u; n<2u*V*L1; ++n, ++X, ++Y) { *Y = *X; }
            X -= 2u*N;
            for (size_t n=0u; n<2u*V*L2; ++n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                X += 2u*K*L2;
                for (size_t b=0u; b<B; ++b, X-=2u*K*L2-2u, Y-=2u*K*L-2u)
                {
                    for (size_t l=0u; l<L1; ++l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
                    X -= 2u*K*L;
                    for (size_t l=0u; l<L2; ++l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
                }
            }
        }
    }
    
    return 0;
}


int cshift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in cshift_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0u) {}
    else if (D==0)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        X += 2u*L2;
        for (size_t l=0u; l<2u*L1; ++l, ++X, ++Y) { *Y = *X; }
        X -= 2u*L;
        for (size_t l=0u; l<2u*L2; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += 2u*L2;
            for (size_t v=0u; v<V; ++v, X+=2u*L)
            {
                for (size_t l=0u; l<2u*L1; ++l, ++X, ++Y) { *Y = *X; }
                X -= 2u*L;
                for (size_t l=0u; l<2u*L2; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            X += 2u*V*L2;
            for (size_t n=0u; n<2u*V*L1; ++n, ++X, ++Y) { *Y = *X; }
            X -= 2u*N;
            for (size_t n=0u; n<2u*V*L2; ++n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                X += 2u*K*L2;
                for (size_t b=0u; b<B; ++b, X-=2u*K*L2-2u, Y-=2u*K*L-2u)
                {
                    for (size_t l=0u; l<L1; ++l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
                    X -= 2u*K*L;
                    for (size_t l=0u; l<L2; ++l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
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
