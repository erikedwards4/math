//Vecs2scalar operation for 2 inputs X1 and X2.
//Kendall rank correlation coefficient for each pair of vectors.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int kendall_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int kendall_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);


int kendall_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in kendall_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in kendall_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float den = (float)(L*(L-1u));
    int s1, s2, ssm = 0;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 1.0f; }
    }
    else if (L==N)
    {
        for (size_t l2=1u; l2<L; ++l2)
        {
            for (size_t l1=0u; l1<l2; ++l1)
            {
                s1 = (X1[l1]>X1[l2]) ? 1 : (X1[l1]<X1[l2]) ? -1 : 0;
                s2 = (X2[l1]>X2[l2]) ? 1 : (X2[l1]<X2[l2]) ? -1 : 0;
                ssm += s1 * s2;
            }
        }
        *Y = 2 * ssm / den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? 0u : L, J2 = (L==N2) ? 0u : L;
            for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, ++Y)
            {
                ssm = 0;
                for (size_t l2=1u; l2<L; ++l2)
                {
                    for (size_t l1=0u; l1<l2; ++l1)
                    {
                        s1 = (X1[l1]>X1[l2]) ? 1 : (X1[l1]<X1[l2]) ? -1 : 0;
                        s2 = (X2[l1]>X2[l2]) ? 1 : (X2[l1]<X2[l2]) ? -1 : 0;
                        ssm += s1 * s2;
                    }
                }
                *Y = 2 * ssm / den;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : B*(L-1u), I2 = (L==N2) ? 0u : B*(L-1u);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1+=J1, X2+=J2, ++Y)
                {
                    ssm = 0;
                    for (size_t l2=1u; l2<L; ++l2)
                    {
                        for (size_t l1=0u; l1<l2; ++l1)
                        {
                            s1 = (X1[l1*K1]>X1[l2*K1]) ? 1 : (X1[l1*K1]<X1[l2*K1]) ? -1 : 0;
                            s2 = (X2[l1*K2]>X2[l2*K2]) ? 1 : (X2[l1*K2]<X2[l2*K2]) ? -1 : 0;
                            ssm += s1 * s2;
                        }
                    }
                    *Y = 2 * ssm / den;
                }
            }
        }
    }

    return 0;
}


int kendall_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in kendall_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in kendall_d: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double den = (double)(L*(L-1u));
    int s1, s2, ssm = 0;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 1.0; }
    }
    else if (L==N)
    {
        for (size_t l2=1u; l2<L; ++l2)
        {
            for (size_t l1=0u; l1<l2; ++l1)
            {
                s1 = (X1[l1]>X1[l2]) ? 1 : (X1[l1]<X1[l2]) ? -1 : 0;
                s2 = (X2[l1]>X2[l2]) ? 1 : (X2[l1]<X2[l2]) ? -1 : 0;
                ssm += s1 * s2;
            }
        }
        *Y = 2 * ssm / den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? 0u : L, J2 = (L==N2) ? 0u : L;
            for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, ++Y)
            {
                ssm = 0;
                for (size_t l2=1u; l2<L; ++l2)
                {
                    for (size_t l1=0u; l1<l2; ++l1)
                    {
                        s1 = (X1[l1]>X1[l2]) ? 1 : (X1[l1]<X1[l2]) ? -1 : 0;
                        s2 = (X2[l1]>X2[l2]) ? 1 : (X2[l1]<X2[l2]) ? -1 : 0;
                        ssm += s1 * s2;
                    }
                }
                *Y = 2 * ssm / den;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : B*(L-1u), I2 = (L==N2) ? 0u : B*(L-1u);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1+=J1, X2+=J2, ++Y)
                {
                    ssm = 0;
                    for (size_t l2=1u; l2<L; ++l2)
                    {
                        for (size_t l1=0u; l1<l2; ++l1)
                        {
                            s1 = (X1[l1*K1]>X1[l2*K1]) ? 1 : (X1[l1*K1]<X1[l2*K1]) ? -1 : 0;
                            s2 = (X2[l1*K2]>X2[l2*K2]) ? 1 : (X2[l1*K2]<X2[l2*K2]) ? -1 : 0;
                            ssm += s1 * s2;
                        }
                    }
                    *Y = 2 * ssm / den;
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
