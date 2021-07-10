//Shifts X by P elements along dim.
//The P edge samples are replaced by 0s.
//This has in-place and not-in-place versions.

//Remarkably, the for loop is faster, even for just the copy 0 part,
//and even for N=200000 (and is much faster for small N)


#include <stdio.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int shift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int shift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int shift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int shift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);

int shift_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int shift_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int shift_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int shift_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);


int shift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (P<0) ? 0 : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0 : (size_t)(-P);

    if (N==0u) {}
    else if (P==0)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        X += L3;
        for (size_t l=0u; l<L1; ++l, ++Y) { *Y = 0.0f; }
        for (size_t l=0u; l<L2; ++l, ++X, ++Y) { *Y = *X; }
        for (size_t l=0u; l<L3; ++l, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            X += L3;
            for (size_t v=0u; v<V; ++v, X+=L1+L3)
            {
                for (size_t l=0u; l<L1; ++l, ++Y) { *Y = 0.0f; }
                for (size_t l=0u; l<L2; ++l, ++X, ++Y) { *Y = *X; }
                for (size_t l=0u; l<L3; ++l, ++Y) { *Y = 0.0f; }
            }
        }
        else if (G==1)
        {
            X += V*L3;
            for (size_t n=0u; n<V*L1; ++n, ++Y) { *Y = 0.0f; }
            for (size_t n=0u; n<V*L2; ++n, ++X, ++Y) { *Y = *X; }
            for (size_t n=0u; n<V*L3; ++n, ++Y) { *Y = 0.0f; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                X += K*L3;
                for (size_t b=0u; b<B; ++b, X-=K*L2-1, Y-=K*L-1)
                {
                    for (size_t l=0u; l<L1; ++l, Y+=K) { *Y = 0.0f; }
                    for (size_t l=0u; l<L2; ++l, X+=K, Y+=K) { *Y = *X; }
                    for (size_t l=0u; l<L3; ++l, Y+=K) { *Y = 0.0f; }
                }
            }
        }
    }
    
    return 0;
}


int shift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (P<0) ? 0 : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0 : (size_t)(-P);

    if (N==0u) {}
    else if (P==0)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        X += L3;
        for (size_t l=0u; l<L1; ++l, ++Y) { *Y = 0.0; }
        for (size_t l=0u; l<L2; ++l, ++X, ++Y) { *Y = *X; }
        for (size_t l=0u; l<L3; ++l, ++Y) { *Y = 0.0; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            X += L3;
            for (size_t v=0u; v<V; ++v, X+=L1+L3)
            {
                for (size_t l=0u; l<L1; ++l, ++Y) { *Y = 0.0; }
                for (size_t l=0u; l<L2; ++l, ++X, ++Y) { *Y = *X; }
                for (size_t l=0u; l<L3; ++l, ++Y) { *Y = 0.0; }
            }
        }
        else if (G==1)
        {
            X += V*L3;
            for (size_t n=0u; n<V*L1; ++n, ++Y) { *Y = 0.0; }
            for (size_t n=0u; n<V*L2; ++n, ++X, ++Y) { *Y = *X; }
            for (size_t n=0u; n<V*L3; ++n, ++Y) { *Y = 0.0; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1), Y+=B*(L-1))
            {
                X += K*L3;
                for (size_t b=0u; b<B; ++b, X-=K*L2-1, Y-=K*L-1)
                {
                    for (size_t l=0u; l<L1; ++l, Y+=K) { *Y = 0.0; }
                    for (size_t l=0u; l<L2; ++l, X+=K, Y+=K) { *Y = *X; }
                    for (size_t l=0u; l<L3; ++l, Y+=K) { *Y = 0.0; }
                }
            }
        }
    }
    
    return 0;
}


int shift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (P<0) ? 0 : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0 : (size_t)(-P);

    if (N==0u) {}
    else if (P==0)
    {
        for (size_t n=0u; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=0u; n<2*N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        X += 2*L3;
        for (size_t l=0u; l<2*L1; ++l, ++Y) { *Y = 0.0f; }
        for (size_t l=0u; l<2*L2; ++l, ++X, ++Y) { *Y = *X; }
        for (size_t l=0u; l<2*L3; ++l, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            X += 2*L3;
            for (size_t v=0u; v<V; ++v, X+=2*(L1+L3))
            {
                for (size_t l=0u; l<2*L1; ++l, ++Y) { *Y = 0.0f; }
                for (size_t l=0u; l<2*L2; ++l, ++X, ++Y) { *Y = *X; }
                for (size_t l=0u; l<2*L3; ++l, ++Y) { *Y = 0.0f; }
            }
        }
        else if (G==1)
        {
            X += 2*V*L3;
            for (size_t n=0u; n<2*V*L1; ++n, ++Y) { *Y = 0.0f; }
            for (size_t n=0u; n<2*V*L2; ++n, ++X, ++Y) { *Y = *X; }
            for (size_t n=0u; n<2*V*L3; ++n, ++Y) { *Y = 0.0f; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1), Y+=2*B*(L-1))
            {
                X += 2*K*L3;
                for (size_t b=0u; b<B; ++b, X-=2*K*L2-2, Y-=2*K*L-2)
                {
                    for (size_t l=0u; l<L1; ++l, Y+=2*K-1) { *Y = 0.0f; *++Y = 0.0f; }
                    for (size_t l=0u; l<L2; ++l, X+=2*K-1, Y+=2*K-1) { *Y = *X; *++Y = *++X; }
                    for (size_t l=0u; l<L3; ++l, Y+=2*K-1) { *Y = 0.0f; *++Y = 0.0f; }
                }
            }
        }
    }
    
    return 0;
}


int shift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (P<0) ? 0 : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0 : (size_t)(-P);

    if (N==0u) {}
    else if (P==0)
    {
        for (size_t n=0u; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=0u; n<2*N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        X += 2*L3;
        for (size_t l=0u; l<2*L1; ++l, ++Y) { *Y = 0.0; }
        for (size_t l=0u; l<2*L2; ++l, ++X, ++Y) { *Y = *X; }
        for (size_t l=0u; l<2*L3; ++l, ++Y) { *Y = 0.0; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            X += 2*L3;
            for (size_t v=0u; v<V; ++v, X+=2*(L1+L3))
            {
                for (size_t l=0u; l<2*L1; ++l, ++Y) { *Y = 0.0; }
                for (size_t l=0u; l<2*L2; ++l, ++X, ++Y) { *Y = *X; }
                for (size_t l=0u; l<2*L3; ++l, ++Y) { *Y = 0.0; }
            }
        }
        else if (G==1)
        {
            X += 2*V*L3;
            for (size_t n=0u; n<2*V*L1; ++n, ++Y) { *Y = 0.0; }
            for (size_t n=0u; n<2*V*L2; ++n, ++X, ++Y) { *Y = *X; }
            for (size_t n=0u; n<2*V*L3; ++n, ++Y) { *Y = 0.0; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1), Y+=2*B*(L-1))
            {
                X += 2*K*L3;
                for (size_t b=0u; b<B; ++b, X-=2*K*L2-2, Y-=2*K*L-2)
                {
                    for (size_t l=0u; l<L1; ++l, Y+=2*K-1) { *Y = 0.0; *++Y = 0.0; }
                    for (size_t l=0u; l<L2; ++l, X+=2*K-1, Y+=2*K-1) { *Y = *X; *++Y = *++X; }
                    for (size_t l=0u; l<L3; ++l, Y+=2*K-1) { *Y = 0.0; *++Y = 0.0; }
                }
            }
        }
    }
    
    return 0;
}


int shift_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (P<0) ? 0 : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0 : (size_t)(-P);

    if (N==0 || P==0) {}
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        if (P<0)
        {
            for (size_t l=0u; l<L2; ++l, ++X) { *X = *(X-P); }
            for (size_t l=0u; l<L3; ++l, ++X) { *X = 0.0f; }
        }
        else
        {
            X += L-1;
            for (size_t l=0u; l<L2; ++l, --X) { *X = *(X-P); }
            for (size_t l=0u; l<L1; ++l, --X) { *X = 0.0f; }
        }
    }
    else if (P<0)
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<L2; ++l, ++X) { *X = *(X-P); }
                for (size_t l=0u; l<L3; ++l, ++X) { *X = 0.0f; }
            }
        }
        else if (G==1)
        {
            for (size_t n=0u; n<V*L2; ++n, ++X) { *X = *(X-P*(int)V); }
            for (size_t n=0u; n<V*L3; ++n, ++X) { *X = 0.0f; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1)
                {
                    for (size_t l=0u; l<L2; ++l, X+=K) { *X = *(X-P*(int)K); }
                    for (size_t l=0u; l<L3; ++l, X+=K) { *X = 0.0f; }
                }
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            X += N-1;
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<L2; ++l, --X) { *X = *(X-P); }
                for (size_t l=0u; l<L1; ++l, --X) { *X = 0.0f; }
            }
        }
        else if (G==1)
        {
            X += N-1;
            for (size_t n=0u; n<V*L2; ++n, --X) { *X = *(X-P*(int)V); }
            for (size_t n=0u; n<V*L1; ++n, --X) { *X = 0.0f; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X+=K+1)
                {
                    X += K*(L-1);
                    for (size_t l=0u; l<L2; ++l, X-=K) { *X = *(X-P*(int)K); }
                    for (size_t l=0u; l<L1; ++l, X-=K) { *X = 0.0f; }
                }
            }
        }
    }
    
    return 0;
}


int shift_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (P<0) ? 0 : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0 : (size_t)(-P);

    if (N==0 || P==0) {}
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        if (P<0)
        {
            for (size_t l=0u; l<L2; ++l, ++X) { *X = *(X-P); }
            for (size_t l=0u; l<L3; ++l, ++X) { *X = 0.0; }
        }
        else
        {
            X += L-1;
            for (size_t l=0u; l<L2; ++l, --X) { *X = *(X-P); }
            for (size_t l=0u; l<L1; ++l, --X) { *X = 0.0; }
        }
    }
    else if (P<0)
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<L2; ++l, ++X) { *X = *(X-P); }
                for (size_t l=0u; l<L3; ++l, ++X) { *X = 0.0; }
            }
        }
        else if (G==1)
        {
            for (size_t n=0u; n<V*L2; ++n, ++X) { *X = *(X-P*(int)V); }
            for (size_t n=0u; n<V*L3; ++n, ++X) { *X = 0.0; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1)
                {
                    for (size_t l=0u; l<L2; ++l, X+=K) { *X = *(X-P*(int)K); }
                    for (size_t l=0u; l<L3; ++l, X+=K) { *X = 0.0; }
                }
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            X += N-1;
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<L2; ++l, --X) { *X = *(X-P); }
                for (size_t l=0u; l<L1; ++l, --X) { *X = 0.0; }
            }
        }
        else if (G==1)
        {
            X += N-1;
            for (size_t n=0u; n<V*L2; ++n, --X) { *X = *(X-P*(int)V); }
            for (size_t n=0u; n<V*L1; ++n, --X) { *X = 0.0; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X+=K+1)
                {
                    X += K*(L-1);
                    for (size_t l=0u; l<L2; ++l, X-=K) { *X = *(X-P*(int)K); }
                    for (size_t l=0u; l<L1; ++l, X-=K) { *X = 0.0; }
                }
            }
        }
    }
    
    return 0;
}


int shift_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_inplace_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (P<0) ? 0 : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0 : (size_t)(-P);

    if (N==0 || P==0) {}
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=0u; n<2*N; ++n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        if (P<0)
        {
            for (size_t l=0u; l<2*L2; ++l, ++X) { *X = *(X-2*P); }
            for (size_t l=0u; l<2*L3; ++l, ++X) { *X = 0.0f; }
        }
        else
        {
            X += 2*L-1;
            for (size_t l=0u; l<2*L2; ++l, --X) { *X = *(X-2*P); }
            for (size_t l=0u; l<2*L1; ++l, --X) { *X = 0.0f; }
        }
    }
    else if (P<0)
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<2*L2; ++l, ++X) { *X = *(X+2*(-P)); }
                for (size_t l=0u; l<2*L3; ++l, ++X) { *X = 0.0f; }
            }
        }
        else if (G==1)
        {
            for (size_t n=0u; n<2*V*L2; ++n, ++X) { *X = *(X-2*P*(int)V); }
            for (size_t n=0u; n<2*V*L3; ++n, ++X) { *X = 0.0f; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2)
                {
                    for (size_t l=0u; l<L2; ++l, X+=2*K) { *X = *(X-2*P*(int)K); *(X+1) = *(X-2*P*(int)K+1); }
                    for (size_t l=0u; l<L3; ++l, X+=2*K) { *X = 0.0f; *(X+1) = 0.0f; }
                }
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            X += 2*N-1;
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<2*L2; ++l, --X) { *X = *(X-2*P); }
                for (size_t l=0u; l<2*L1; ++l, --X) { *X = 0.0f; }
            }
        }
        else if (G==1)
        {
            X += 2*N-1;
            for (size_t n=0u; n<2*V*L2; ++n, --X) { *X = *(X-2*P*(int)V); }
            for (size_t n=0u; n<2*V*L1; ++n, --X) { *X = 0.0f; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X+=2*K+2)
                {
                    X += 2*K*(L-1);
                    for (size_t l=0u; l<L2; ++l, X-=2*K) { *X = *(X-2*P*(int)K); *(X+1) = *(X-2*P*(int)K+1); }
                    for (size_t l=0u; l<L1; ++l, X-=2*K) { *X = 0.0f; *(X+1) = 0.0f; }
                }
            }
        }
    }
    
    return 0;
}


int shift_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_inplace_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (P<0) ? 0 : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0 : (size_t)(-P);

    if (N==0 || P==0) {}
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=0u; n<2*N; ++n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        if (P<0)
        {
            for (size_t l=0u; l<2*L2; ++l, ++X) { *X = *(X-2*P); }
            for (size_t l=0u; l<2*L3; ++l, ++X) { *X = 0.0; }
        }
        else
        {
            X += 2*L-1;
            for (size_t l=0u; l<2*L2; ++l, --X) { *X = *(X-2*P); }
            for (size_t l=0u; l<2*L1; ++l, --X) { *X = 0.0; }
        }
    }
    else if (P<0)
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<2*L2; ++l, ++X) { *X = *(X+2*(-P)); }
                for (size_t l=0u; l<2*L3; ++l, ++X) { *X = 0.0; }
            }
        }
        else if (G==1)
        {
            for (size_t n=0u; n<2*V*L2; ++n, ++X) { *X = *(X-2*P*(int)V); }
            for (size_t n=0u; n<2*V*L3; ++n, ++X) { *X = 0.0; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2)
                {
                    for (size_t l=0u; l<L2; ++l, X+=2*K) { *X = *(X-2*P*(int)K); *(X+1) = *(X-2*P*(int)K+1); }
                    for (size_t l=0u; l<L3; ++l, X+=2*K) { *X = 0.0; *(X+1) = 0.0; }
                }
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            X += 2*N-1;
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<2*L2; ++l, --X) { *X = *(X-2*P); }
                for (size_t l=0u; l<2*L1; ++l, --X) { *X = 0.0; }
            }
        }
        else if (G==1)
        {
            X += 2*N-1;
            for (size_t n=0u; n<2*V*L2; ++n, --X) { *X = *(X-2*P*(int)V); }
            for (size_t n=0u; n<2*V*L1; ++n, --X) { *X = 0.0; }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X+=2*K+2)
                {
                    X += 2*K*(L-1);
                    for (size_t l=0u; l<L2; ++l, X-=2*K) { *X = *(X-2*P*(int)K); *(X+1) = *(X-2*P*(int)K+1); }
                    for (size_t l=0u; l<L1; ++l, X-=2*K) { *X = 0.0; *(X+1) = 0.0; }
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
