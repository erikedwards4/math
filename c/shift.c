//Shifts X by P elements along dim.
//The P edge samples are replaced by 0s.
//This has in-place and not-in-place versions.

//Remarkably, the for loop is faster, even for just the copy 0 part,
//and even for N=200000 (and is much faster for small N)


#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int shift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in shift_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (P<0) ? 0u : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0u : (size_t)(-P);

    if (N==0u) {}
    else if (P==0)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        X += L3;
        for (size_t l=L1; l>0u; --l, ++Y) { *Y = 0.0f; }
        for (size_t l=L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
        for (size_t l=L3; l>0u; --l, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += L3;
            for (size_t v=V; v>0u; --v, X+=L1+L3)
            {
                for (size_t l=L1; l>0u; --l, ++Y) { *Y = 0.0f; }
                for (size_t l=L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
                for (size_t l=L3; l>0u; --l, ++Y) { *Y = 0.0f; }
            }
        }
        else if (G==1u)
        {
            X += V*L3;
            for (size_t n=V*L1; n>0u; --n, ++Y) { *Y = 0.0f; }
            for (size_t n=V*L2; n>0u; --n, ++X, ++Y) { *Y = *X; }
            for (size_t n=V*L3; n>0u; --n, ++Y) { *Y = 0.0f; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                X += K*L3;
                for (size_t b=B; b>0u; --b, X-=K*L2-1u, Y-=K*L-1u)
                {
                    for (size_t l=L1; l>0u; --l, Y+=K) { *Y = 0.0f; }
                    for (size_t l=L2; l>0u; --l, X+=K, Y+=K) { *Y = *X; }
                    for (size_t l=L3; l>0u; --l, Y+=K) { *Y = 0.0f; }
                }
            }
        }
    }
    
    return 0;
}


int shift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in shift_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (P<0) ? 0u : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0u : (size_t)(-P);

    if (N==0u) {}
    else if (P==0)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        X += L3;
        for (size_t l=L1; l>0u; --l, ++Y) { *Y = 0.0; }
        for (size_t l=L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
        for (size_t l=L3; l>0u; --l, ++Y) { *Y = 0.0; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += L3;
            for (size_t v=V; v>0u; --v, X+=L1+L3)
            {
                for (size_t l=L1; l>0u; --l, ++Y) { *Y = 0.0; }
                for (size_t l=L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
                for (size_t l=L3; l>0u; --l, ++Y) { *Y = 0.0; }
            }
        }
        else if (G==1u)
        {
            X += V*L3;
            for (size_t n=V*L1; n>0u; --n, ++Y) { *Y = 0.0; }
            for (size_t n=V*L2; n>0u; --n, ++X, ++Y) { *Y = *X; }
            for (size_t n=V*L3; n>0u; --n, ++Y) { *Y = 0.0; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                X += K*L3;
                for (size_t b=B; b>0u; --b, X-=K*L2-1u, Y-=K*L-1u)
                {
                    for (size_t l=L1; l>0u; --l, Y+=K) { *Y = 0.0; }
                    for (size_t l=L2; l>0u; --l, X+=K, Y+=K) { *Y = *X; }
                    for (size_t l=L3; l>0u; --l, Y+=K) { *Y = 0.0; }
                }
            }
        }
    }
    
    return 0;
}


int shift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in shift_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (P<0) ? 0u : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0u : (size_t)(-P);

    if (N==0u) {}
    else if (P==0)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=2u*N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        X += 2u*L3;
        for (size_t l=2u*L1; l>0u; --l, ++Y) { *Y = 0.0f; }
        for (size_t l=2u*L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
        for (size_t l=2u*L3; l>0u; --l, ++Y) { *Y = 0.0f; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += 2u*L3;
            for (size_t v=V; v>0u; --v, X+=2u*(L1+L3))
            {
                for (size_t l=2u*L1; l>0u; --l, ++Y) { *Y = 0.0f; }
                for (size_t l=2u*L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
                for (size_t l=2u*L3; l>0u; --l, ++Y) { *Y = 0.0f; }
            }
        }
        else if (G==1u)
        {
            X += 2u*V*L3;
            for (size_t n=2u*V*L1; n>0u; --n, ++Y) { *Y = 0.0f; }
            for (size_t n=2u*V*L2; n>0u; --n, ++X, ++Y) { *Y = *X; }
            for (size_t n=2u*V*L3; n>0u; --n, ++Y) { *Y = 0.0f; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                X += 2u*K*L3;
                for (size_t b=B; b>0u; --b, X-=2u*K*L2-2u, Y-=2u*K*L-2u)
                {
                    for (size_t l=L1; l>0u; --l, Y+=2u*K-1u) { *Y = 0.0f; *++Y = 0.0f; }
                    for (size_t l=L2; l>0u; --l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
                    for (size_t l=L3; l>0u; --l, Y+=2u*K-1u) { *Y = 0.0f; *++Y = 0.0f; }
                }
            }
        }
    }
    
    return 0;
}


int shift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in shift_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (P<0) ? 0u : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0u : (size_t)(-P);

    if (N==0u) {}
    else if (P==0)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=2u*N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        X += 2u*L3;
        for (size_t l=2u*L1; l>0u; --l, ++Y) { *Y = 0.0; }
        for (size_t l=2u*L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
        for (size_t l=2u*L3; l>0u; --l, ++Y) { *Y = 0.0; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += 2u*L3;
            for (size_t v=V; v>0u; --v, X+=2u*(L1+L3))
            {
                for (size_t l=2u*L1; l>0u; --l, ++Y) { *Y = 0.0; }
                for (size_t l=2u*L2; l>0u; --l, ++X, ++Y) { *Y = *X; }
                for (size_t l=2u*L3; l>0u; --l, ++Y) { *Y = 0.0; }
            }
        }
        else if (G==1u)
        {
            X += 2u*V*L3;
            for (size_t n=2u*V*L1; n>0u; --n, ++Y) { *Y = 0.0; }
            for (size_t n=2u*V*L2; n>0u; --n, ++X, ++Y) { *Y = *X; }
            for (size_t n=2u*V*L3; n>0u; --n, ++Y) { *Y = 0.0; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                X += 2u*K*L3;
                for (size_t b=B; b>0u; --b, X-=2u*K*L2-2u, Y-=2u*K*L-2u)
                {
                    for (size_t l=L1; l>0u; --l, Y+=2u*K-1u) { *Y = 0.0; *++Y = 0.0; }
                    for (size_t l=L2; l>0u; --l, X+=2u*K-1u, Y+=2u*K-1u) { *Y = *X; *++Y = *++X; }
                    for (size_t l=L3; l>0u; --l, Y+=2u*K-1u) { *Y = 0.0; *++Y = 0.0; }
                }
            }
        }
    }
    
    return 0;
}


int shift_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in shift_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (P<0) ? 0 : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0 : (size_t)(-P);

    if (N==0u || P==0) {}
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        if (P<0)
        {
            for (size_t l=L2; l>0u; --l, ++X) { *X = *(X-P); }
            for (size_t l=L3; l>0u; --l, ++X) { *X = 0.0f; }
        }
        else
        {
            X += L-1u;
            for (size_t l=L2; l>0u; --l, --X) { *X = *(X-P); }
            for (size_t l=L1; l>0u; --l, --X) { *X = 0.0f; }
        }
    }
    else if (P<0)
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=L2; l>0u; --l, ++X) { *X = *(X-P); }
                for (size_t l=L3; l>0u; --l, ++X) { *X = 0.0f; }
            }
        }
        else if (G==1u)
        {
            for (size_t n=V*L2; n>0u; --n, ++X) { *X = *(X-P*(int)V); }
            for (size_t n=V*L3; n>0u; --n, ++X) { *X = 0.0f; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u)
                {
                    for (size_t l=L2; l>0u; --l, X+=K) { *X = *(X-P*(int)K); }
                    for (size_t l=L3; l>0u; --l, X+=K) { *X = 0.0f; }
                }
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += N-1u;
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=L2; l>0u; --l, --X) { *X = *(X-P); }
                for (size_t l=L1; l>0u; --l, --X) { *X = 0.0f; }
            }
        }
        else if (G==1u)
        {
            X += N-1u;
            for (size_t n=V*L2; n>0u; --n, --X) { *X = *(X-P*(int)V); }
            for (size_t n=V*L1; n>0u; --n, --X) { *X = 0.0f; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X+=K+1u)
                {
                    X += K*(L-1u);
                    for (size_t l=L2; l>0u; --l, X-=K) { *X = *(X-P*(int)K); }
                    for (size_t l=L1; l>0u; --l, X-=K) { *X = 0.0f; }
                }
            }
        }
    }
    
    return 0;
}


int shift_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in shift_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (P<0) ? 0 : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0 : (size_t)(-P);

    if (N==0u || P==0) {}
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=N; n>0u; --n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        if (P<0)
        {
            for (size_t l=L2; l>0u; --l, ++X) { *X = *(X-P); }
            for (size_t l=L3; l>0u; --l, ++X) { *X = 0.0; }
        }
        else
        {
            X += L-1u;
            for (size_t l=L2; l>0u; --l, --X) { *X = *(X-P); }
            for (size_t l=L1; l>0u; --l, --X) { *X = 0.0; }
        }
    }
    else if (P<0)
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=L2; l>0u; --l, ++X) { *X = *(X-P); }
                for (size_t l=L3; l>0u; --l, ++X) { *X = 0.0; }
            }
        }
        else if (G==1u)
        {
            for (size_t n=V*L2; n>0u; --n, ++X) { *X = *(X-P*(int)V); }
            for (size_t n=V*L3; n>0u; --n, ++X) { *X = 0.0; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u)
                {
                    for (size_t l=L2; l>0u; --l, X+=K) { *X = *(X-P*(int)K); }
                    for (size_t l=L3; l>0u; --l, X+=K) { *X = 0.0; }
                }
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += N-1u;
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=L2; l>0u; --l, --X) { *X = *(X-P); }
                for (size_t l=L1; l>0u; --l, --X) { *X = 0.0; }
            }
        }
        else if (G==1u)
        {
            X += N-1u;
            for (size_t n=V*L2; n>0u; --n, --X) { *X = *(X-P*(int)V); }
            for (size_t n=V*L1; n>0u; --n, --X) { *X = 0.0; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X+=K+1u)
                {
                    X += K*(L-1u);
                    for (size_t l=L2; l>0u; --l, X-=K) { *X = *(X-P*(int)K); }
                    for (size_t l=L1; l>0u; --l, X-=K) { *X = 0.0; }
                }
            }
        }
    }
    
    return 0;
}


int shift_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in shift_inplace_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (P<0) ? 0 : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0 : (size_t)(-P);

    if (N==0u || P==0) {}
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=2u*N; n>0u; --n, ++X) { *X = 0.0f; }
    }
    else if (L==N)
    {
        if (P<0)
        {
            for (size_t l=2u*L2; l>0u; --l, ++X) { *X = *(X-2*P); }
            for (size_t l=2u*L3; l>0u; --l, ++X) { *X = 0.0f; }
        }
        else
        {
            X += 2u*L-1u;
            for (size_t l=2u*L2; l>0u; --l, --X) { *X = *(X-2*P); }
            for (size_t l=2u*L1; l>0u; --l, --X) { *X = 0.0f; }
        }
    }
    else if (P<0)
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=2u*L2; l>0u; --l, ++X) { *X = *(X+2*(-P)); }
                for (size_t l=2u*L3; l>0u; --l, ++X) { *X = 0.0f; }
            }
        }
        else if (G==1u)
        {
            for (size_t n=2u*V*L2; n>0u; --n, ++X) { *X = *(X-2*P*(int)V); }
            for (size_t n=2u*V*L3; n>0u; --n, ++X) { *X = 0.0f; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u)
                {
                    for (size_t l=L2; l>0u; --l, X+=2u*K) { *X = *(X-2*P*(int)K); *(X+1) = *(X-2*P*(int)K+1u); }
                    for (size_t l=L3; l>0u; --l, X+=2u*K) { *X = 0.0f; *(X+1) = 0.0f; }
                }
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += 2u*N-1u;
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=2u*L2; l>0u; --l, --X) { *X = *(X-2*P); }
                for (size_t l=2u*L1; l>0u; --l, --X) { *X = 0.0f; }
            }
        }
        else if (G==1u)
        {
            X += 2u*N-1u;
            for (size_t n=2u*V*L2; n>0u; --n, --X) { *X = *(X-2*P*(int)V); }
            for (size_t n=2u*V*L1; n>0u; --n, --X) { *X = 0.0f; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X+=2u*K+2u)
                {
                    X += 2u*K*(L-1u);
                    for (size_t l=L2; l>0u; --l, X-=2u*K) { *X = *(X-2*P*(int)K); *(X+1) = *(X-2*P*(int)K+1u); }
                    for (size_t l=L1; l>0u; --l, X-=2u*K) { *X = 0.0f; *(X+1) = 0.0f; }
                }
            }
        }
    }
    
    return 0;
}


int shift_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int P)
{
    if (dim>3u) { fprintf(stderr,"error in shift_inplace_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (P<0) ? 0 : (size_t)P;
    const size_t L2 = (P<0) ? L-(size_t)(-P) : L-(size_t)P;
    const size_t L3 = (P>0) ? 0 : (size_t)(-P);

    if (N==0u || P==0) {}
    else if (P<=-(int)L || P>=(int)L)
    {
        for (size_t n=2u*N; n>0u; --n, ++X) { *X = 0.0; }
    }
    else if (L==N)
    {
        if (P<0)
        {
            for (size_t l=2u*L2; l>0u; --l, ++X) { *X = *(X-2*P); }
            for (size_t l=2u*L3; l>0u; --l, ++X) { *X = 0.0; }
        }
        else
        {
            X += 2u*L-1u;
            for (size_t l=2u*L2; l>0u; --l, --X) { *X = *(X-2*P); }
            for (size_t l=2u*L1; l>0u; --l, --X) { *X = 0.0; }
        }
    }
    else if (P<0)
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=2u*L2; l>0u; --l, ++X) { *X = *(X+2*(-P)); }
                for (size_t l=2u*L3; l>0u; --l, ++X) { *X = 0.0; }
            }
        }
        else if (G==1u)
        {
            for (size_t n=2u*V*L2; n>0u; --n, ++X) { *X = *(X-2*P*(int)V); }
            for (size_t n=2u*V*L3; n>0u; --n, ++X) { *X = 0.0; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u)
                {
                    for (size_t l=L2; l>0u; --l, X+=2u*K) { *X = *(X-2*P*(int)K); *(X+1) = *(X-2*P*(int)K+1u); }
                    for (size_t l=L3; l>0u; --l, X+=2u*K) { *X = 0.0; *(X+1) = 0.0; }
                }
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            X += 2u*N-1u;
            for (size_t v=V; v>0u; --v)
            {
                for (size_t l=2u*L2; l>0u; --l, --X) { *X = *(X-2*P); }
                for (size_t l=2u*L1; l>0u; --l, --X) { *X = 0.0; }
            }
        }
        else if (G==1u)
        {
            X += 2u*N-1u;
            for (size_t n=2u*V*L2; n>0u; --n, --X) { *X = *(X-2*P*(int)V); }
            for (size_t n=2u*V*L1; n>0u; --n, --X) { *X = 0.0; }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X+=2u*K+2u)
                {
                    X += 2u*K*(L-1u);
                    for (size_t l=L2; l>0u; --l, X-=2u*K) { *X = *(X-2*P*(int)K); *(X+1) = *(X-2*P*(int)K+1u); }
                    for (size_t l=L1; l>0u; --l, X-=2u*K) { *X = 0.0; *(X+1) = 0.0; }
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
