//Vecs2scalar operation for 2 inputs X1 and X2.
//Cosine distance for each pair of vectors.
//This is 1 - cosine similarity, so is in [0 2].
//Although not a proper distance metric, it is often used.

//Broadcasting is only meant for 1D vectors into larger tensors,
//so if using broadcasting, either X1 or X2 must be a vector,
//and the shapes must be equal along the appropiate axis.

#include <stdio.h>
#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int dist_cos_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dist_cos_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist_cos_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    float sm11, sm22, sm12;

    if (N==0u) {}
    else if (L==N)
    {
        sm11 = sm22 = sm12 = 0.0f;
        for (size_t l=L; l>0u; --l, ++X1, ++X2)
        {
            sm11 += *X1**X1; sm22 += *X2**X2; sm12 += *X1**X2;
        }
        *Y = 1.0f - sm12/sqrtf(sm11*sm22);
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
                sm11 = sm22 = sm12 = 0.0f;
                for (size_t l=L; l>0u; --l, ++X1, ++X2)
                {
                    sm11 += *X1**X1; sm22 += *X2**X2; sm12 += *X1**X2;
                }
                *Y = 1.0f - sm12/sqrtf(sm11*sm22);
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
                    sm11 = sm22 = sm12 = 0.0f;
                    for (size_t l=L; l>0u; --l, X1+=K1, X2+=K2)
                    {
                        sm11 += *X1**X1; sm22 += *X2**X2; sm12 += *X1**X2;
                    }
                    *Y = 1.0f - sm12/sqrtf(sm11*sm22);
                }
            }
        }
    }

    return 0;
}


int dist_cos_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dist_cos_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist_cos_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    double sm11, sm22, sm12;

    if (N==0u) {}
    else if (L==N)
    {
        sm11 = sm22 = sm12 = 0.0;
        for (size_t l=L; l>0u; --l, ++X1, ++X2)
        {
            sm11 += *X1**X1; sm22 += *X2**X2; sm12 += *X1**X2;
        }
        *Y = 1.0 - sm12/sqrt(sm11*sm22);
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
                sm11 = sm22 = sm12 = 0.0;
                for (size_t l=L; l>0u; --l, ++X1, ++X2)
                {
                    sm11 += *X1**X1; sm22 += *X2**X2; sm12 += *X1**X2;
                }
                *Y = 1.0 - sm12/sqrt(sm11*sm22);
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
                    sm11 = sm22 = sm12 = 0.0;
                    for (size_t l=L; l>0u; --l, X1+=K1, X2+=K2)
                    {
                        sm11 += *X1**X1; sm22 += *X2**X2; sm12 += *X1**X2;
                    }
                    *Y = 1.0 - sm12/sqrt(sm11*sm22);
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
