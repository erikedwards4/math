//Circular shift of X by P elements along dim.
//This has no in-place version (no faster version that I can find).

#include <stdio.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cshift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int cshift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int cshift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int cshift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);


int cshift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in cshift_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0) {}
    else if (D==0)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        X += L2;
        for (size_t l=0; l<L1; ++l, ++X, ++Y) { *Y = *X; }
        X -= L;
        for (size_t l=0; l<L2; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            X += L2;
            for (size_t v=0; v<V; ++v, X+=L)
            {
                for (size_t l=0; l<L1; ++l, ++X, ++Y) { *Y = *X; }
                X -= L;
                for (size_t l=0; l<L2; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1)
        {
            X += V * L2;
            for (size_t n=0; n<V*L1; ++n, ++X, ++Y) { *Y = *X; }
            X -= N;
            for (size_t n=0; n<V*L2; ++n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                X += K * L2;
                for (size_t b=0; b<B; ++b, X-=K*L2-1, Y-=K*L-1)
                {
                    for (size_t l=0; l<L1; ++l, X+=K, Y+=K) { *Y = *X; }
                    X -= K*L;
                    for (size_t l=0; l<L2; ++l, X+=K, Y+=K) { *Y = *X; }
                }
            }
        }
    }
    
    return 0;
}


int cshift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in cshift_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0) {}
    else if (D==0)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        X += L2;
        for (size_t l=0; l<L1; ++l, ++X, ++Y) { *Y = *X; }
        X -= L;
        for (size_t l=0; l<L2; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            X += L2;
            for (size_t v=0; v<V; ++v, X+=L)
            {
                for (size_t l=0; l<L1; ++l, ++X, ++Y) { *Y = *X; }
                X -= L;
                for (size_t l=0; l<L2; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1)
        {
            X += V * L2;
            for (size_t n=0; n<V*L1; ++n, ++X, ++Y) { *Y = *X; }
            X -= N;
            for (size_t n=0; n<V*L2; ++n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                X += K * L2;
                for (size_t b=0; b<B; ++b, X-=K*L2-1, Y-=K*L-1)
                {
                    for (size_t l=0; l<L1; ++l, X+=K, Y+=K) { *Y = *X; }
                    X -= K*L;
                    for (size_t l=0; l<L2; ++l, X+=K, Y+=K) { *Y = *X; }
                }
            }
        }
    }
    
    return 0;
}


int cshift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in cshift_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0) {}
    else if (D==0)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        X += 2*L2;
        for (size_t l=0; l<2*L1; ++l, ++X, ++Y) { *Y = *X; }
        X -= 2*L;
        for (size_t l=0; l<2*L2; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            X += 2*L2;
            for (size_t v=0; v<V; ++v, X+=2*L)
            {
                for (size_t l=0; l<2*L1; ++l, ++X, ++Y) { *Y = *X; }
                X -= 2*L;
                for (size_t l=0; l<2*L2; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1)
        {
            X += 2*V*L2;
            for (size_t n=0; n<2*V*L1; ++n, ++X, ++Y) { *Y = *X; }
            X -= 2*N;
            for (size_t n=0; n<2*V*L2; ++n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1), Y+=2*B*(L-1))
            {
                X += 2*K*L2;
                for (size_t b=0; b<B; ++b, X-=2*K*L2-2, Y-=2*K*L-2)
                {
                    for (size_t l=0; l<L1; ++l, X+=2*K-1, Y+=2*K-1) { *Y = *X; *++Y = *++X; }
                    X -= 2*K*L;
                    for (size_t l=0; l<L2; ++l, X+=2*K-1, Y+=2*K-1) { *Y = *X; *++Y = *++X; }
                }
            }
        }
    }
    
    return 0;
}


int cshift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in cshift_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    int D = P % (int)L;
    const size_t L1 = (D<0) ? (size_t)((int)L+D) : (size_t)D;
    const size_t L2 = (D<0) ? (size_t)(-D) : L-(size_t)D;

    if (N==0) {}
    else if (D==0)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        X += 2*L2;
        for (size_t l=0; l<2*L1; ++l, ++X, ++Y) { *Y = *X; }
        X -= 2*L;
        for (size_t l=0; l<2*L2; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            X += 2*L2;
            for (size_t v=0; v<V; ++v, X+=2*L)
            {
                for (size_t l=0; l<2*L1; ++l, ++X, ++Y) { *Y = *X; }
                X -= 2*L;
                for (size_t l=0; l<2*L2; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else if (G==1)
        {
            X += 2*V*L2;
            for (size_t n=0; n<2*V*L1; ++n, ++X, ++Y) { *Y = *X; }
            X -= 2*N;
            for (size_t n=0; n<2*V*L2; ++n, ++X, ++Y) { *Y = *X; }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1), Y+=2*B*(L-1))
            {
                X += 2*K*L2;
                for (size_t b=0; b<B; ++b, X-=2*K*L2-2, Y-=2*K*L-2)
                {
                    for (size_t l=0; l<L1; ++l, X+=2*K-1, Y+=2*K-1) { *Y = *X; *++Y = *++X; }
                    X -= 2*K*L;
                    for (size_t l=0; l<L2; ++l, X+=2*K-1, Y+=2*K-1) { *Y = *X; *++Y = *++X; }
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
