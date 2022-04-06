//Vecs2scalar operation for 2 inputs X1 and X2.
//Gets L0 distance between pairs of vecs in X1, X2.

//This is the 'Hamming distance',
//which is just the count of non-zero values in |X1-X2|,
//which is just the count of X1!=X2.
//However, for numerical purposes, I use |X1-X2|>thresh.

//Broadcasting is only meant for 1D vectors into larger tensors,
//so if using broadcasting, either X1 or X2 must be a vector,
//and the shapes must be equal along the appropiate axis.

#include <stdio.h>
#include <float.h>
#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int dist0_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dist0_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist0_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    float thresh = 2.0f * FLT_EPSILON;

    if (N==0u) {}
    else if (L==N)
    {
        int cnt = 0;
        for (size_t l=L; l>0u; --l, ++X1, ++X2)
        {
            cnt += (fabsf(*X1-*X2)>thresh);
        }
        *Y = (float)cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? L : 0u, J2 = (L==N2) ? L : 0u;
            for (size_t v=V; v>0u; --v, X1-=J1, X2-=J2, ++Y)
            {
                int cnt = 0;
                for (size_t l=L; l>0u; --l, ++X1, ++X2)
                {
                    //cnt += (*X1!=*X2);
                    cnt += (fabsf(*X1-*X2)>thresh);
                }
                *Y = (float)cnt;
            }
        }
        else if (G==1u)
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            for (size_t v=V; v>0u; --v, X1+=J1, X2+=J2, ++Y)
            {
                //*Y = (float)(*X1!=*X2);
                *Y = (float)(fabsf(*X1-*X2)>thresh);
            }
            X1 += 1u-J1; X2 += 1u-J2; Y -= V;
            for (size_t l=L; l>1u; --l, X1+=1u-J1, X2+=1u-J2, Y-=V)
            {
                for (size_t v=V; v>0u; --v, X1+=J1, X2+=J2, ++Y)
                {
                    //*Y += (float)(*X1!=*X2);
                    *Y += (float)(fabsf(*X1-*X2)>thresh);
                }
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : B*(L-1u), I2 = (L==N2) ? 0u : B*(L-1u);
            for (size_t g=G; g>0u; --g, X1+=I1, X2+=I2)
            {
                for (size_t b=B; b>0u; --b, X1-=L*K1-J1, X2-=L*K2-J2)
                {
                    int cnt = 0;
                    for (size_t l=L; l>0u; --l, X1+=K1, X2+=K2)
                    {
                        //cnt += (*X1!=*X2);
                        cnt += (fabsf(*X1-*X2)>thresh);
                    }
                    *Y = (float)cnt;
                }
            }
        }
    }

    return 0;
}


int dist0_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dist0_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist0_d: vectors in X1 and X2 must have the same length\n"); return 1; }
    double thresh = 2.0 * DBL_EPSILON;

    if (N==0u) {}
    else if (L==N)
    {
        int cnt = 0;
        for (size_t l=L; l>0u; --l, ++X1, ++X2)
        {
            cnt += (fabs(*X1-*X2)>thresh);
        }
        *Y = (double)cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? L : 0u, J2 = (L==N2) ? L : 0u;
            for (size_t v=V; v>0u; --v, X1-=J1, X2-=J2, ++Y)
            {
                int cnt = 0;
                for (size_t l=L; l>0u; --l, ++X1, ++X2)
                {
                    cnt += (fabs(*X1-*X2)>thresh);
                }
                *Y = (double)cnt;
            }
        }
        else if (G==1u)
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            for (size_t v=V; v>0u; --v, X1+=J1, X2+=J2, ++Y)
            {
                *Y = (double)(fabs(*X1-*X2)>thresh);
            }
            X1 += 1u-J1; X2 += 1u-J2; Y -= V;
            for (size_t l=L; l>1u; --l, X1+=1u-J1, X2+=1u-J2, Y-=V)
            {
                for (size_t v=V; v>0u; --v, X1+=J1, X2+=J2, ++Y)
                {
                    *Y += (double)(fabs(*X1-*X2)>thresh);
                }
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : B*(L-1u), I2 = (L==N2) ? 0u : B*(L-1u);
            for (size_t g=G; g>0u; --g, X1+=I1, X2+=I2)
            {
                for (size_t b=B; b>0u; --b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    int cnt = 0;
                    for (size_t l=L; l>0u; --l, X1+=K1, X2+=K2)
                    {
                        cnt += (fabs(*X1-*X2)>thresh);
                    }
                    *Y = (double)cnt;
                }
            }
        }
    }

    return 0;
}


int dist0_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dist0_c: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist0_c: vectors in X1 and X2 must have the same length\n"); return 1; }
    float thresh = 2.0f * FLT_EPSILON;

    if (N==0u) {}
    else if (L==N)
    {
        int cnt = 0;
        for (size_t l=L; l>0u; --l, ++X1, ++X2)
        {
            cnt += (fabsf(*X1++-*X2++)>thresh || fabsf(*X1-*X2)>thresh);
        }
        *Y = (float)cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? 2u*L : 0u, J2 = (L==N2) ? 2u*L : 0u;
            for (size_t v=V; v>0u; --v, X1-=J1, X2-=J2, ++Y)
            {
                int cnt = 0;
                for (size_t l=L; l>0u; --l, ++X1, ++X2)
                {
                    cnt += (fabsf(*X1++-*X2++)>thresh || fabsf(*X1-*X2)>thresh);
                }
                *Y = (float)cnt;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 2u, J2 = (L==N2) ? 0u : 2u;
            const size_t K1 = (L==N1) ? 2u : 2u*K, K2 = (L==N2) ? 2u : 2u*K;
            const size_t I1 = (L==N1) ? 0u : 2u*B*(L-1u), I2 = (L==N2) ? 0u : 2u*B*(L-1u);
            for (size_t g=G; g>0u; --g, X1+=I1, X2+=I2)
            {
                for (size_t b=B; b>0u; --b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    int cnt = 0;
                    for (size_t l=L; l>0u; --l, X1+=K1-1u, X2+=K2-1u)
                    {
                        cnt += (fabsf(*X1++-*X2++)>thresh || fabsf(*X1-*X2)>thresh);
                    }
                    *Y = (float)cnt;
                }
            }
        }
    }

    return 0;
}


int dist0_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dist0_z: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist0_z: vectors in X1 and X2 must have the same length\n"); return 1; }
    double thresh = 2.0 * DBL_EPSILON;

    if (N==0u) {}
    else if (L==N)
    {
        int cnt = 0;
        for (size_t l=L; l>0u; --l, ++X1, ++X2)
        {
            cnt += (fabs(*X1++-*X2++)>thresh || fabs(*X1-*X2)>thresh);
        }
        *Y = (double)cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? 2u*L : 0u, J2 = (L==N2) ? 2u*L : 0u;
            for (size_t v=V; v>0u; --v, X1-=J1, X2-=J2, ++Y)
            {
                int cnt = 0;
                for (size_t l=L; l>0u; --l, ++X1, ++X2)
                {
                    cnt += (fabs(*X1++-*X2++)>thresh || fabs(*X1-*X2)>thresh);
                }
                *Y = (double)cnt;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 2u, J2 = (L==N2) ? 0u : 2u;
            const size_t K1 = (L==N1) ? 2u : 2u*K, K2 = (L==N2) ? 2u : 2u*K;
            const size_t I1 = (L==N1) ? 0u : 2u*B*(L-1u), I2 = (L==N2) ? 0u : 2u*B*(L-1u);
            for (size_t g=G; g>0u; --g, X1+=I1, X2+=I2)
            {
                for (size_t b=B; b>0u; --b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    int cnt = 0;
                    for (size_t l=L; l>0u; --l, X1+=K1-1u, X2+=K2-1u)
                    {
                        cnt += (fabs(*X1++-*X2++)>thresh || fabs(*X1-*X2)>thresh);
                    }
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
