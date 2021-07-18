//Vecs2scalar operation for 2 inputs X1 and X2.
//Cokurtosis for each pair of vectors.
//This is essentially the corr of squared values.

//For consistency with Wikipedia and the univariate kurtosis function,
//the final result is multiplied by L (so no longer in [0 1]).
//Thus, cokurtosis(X,X) = kurtosis(X), using the biased version of kurtosis.

//This is a measure of co-occurrence of extreme positive and negative deviations.
//It is invariant to scale: cokurtosis(X1,X2) = cokurtosis(a+b*X1,c+d*X2).

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cokurtosis_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int cokurtosis_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);


int cokurtosis_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cokurtosis_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in cokurtosis_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float den = 1.0f / (float)L;
    float xx1, xx2, mn1, mn2, ss1, ss2, ss12;

    if (N==0u) {}
    else if (L==N)
    {
        mn1 = mn2 = ss1 = ss2 = ss12 = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2) { mn1 += *X1; mn2 += *X2; }
        mn1 *= den; mn2 *= den;
        X1 -= L; X2 -= L;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            xx1 = *X1 - mn1; xx1 *= xx1;
            xx2 = *X2 - mn2; xx2 *= xx2;
            ss1 += xx1; ss2 += xx2; ss12 += xx1*xx2;
        }
        *Y = (float)L * ss12 / (ss1*ss2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? L : 0u, J2 = (L==N2) ? L : 0u;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                mn1 = mn2 = ss1 = ss2 = ss12 = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2) { mn1 += *X1; mn2 += *X2; }
                mn1 *= den; mn2 *= den;
                X1 -= L; X2 -= L;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    xx1 = *X1 - mn1; xx1 *= xx1;
                    xx2 = *X2 - mn2; xx2 *= xx2;
                    ss1 += xx1; ss2 += xx2; ss12 += xx1*xx2;
                }
                *Y = (float)L * ss12 / (ss1*ss2);
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : B*(L-1u), I2 = (L==N2) ? 0u : B*(L-1u);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    mn1 = mn2 = ss1 = ss2 = ss12 = 0.0f;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2) { mn1 += *X1; mn2 += *X2; }
                    mn1 *= den; mn2 *= den;
                    X1 -= L*K1; X2 -= L*K2;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        xx1 = *X1 - mn1; xx1 *= xx1;
                        xx2 = *X2 - mn2; xx2 *= xx2;
                        ss1 += xx1; ss2 += xx2; ss12 += xx1*xx2;
                    }
                    *Y = (float)L * ss12 / (ss1*ss2);
                }
            }
        }
    }

    return 0;
}


int cokurtosis_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in cokurtosis_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in cokurtosis_d: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double den = 1.0 / (double)L;
    double xx1, xx2, mn1, mn2, ss1, ss2, ss12;

    if (N==0u) {}
    else if (L==N)
    {
        mn1 = mn2 = ss1 = ss2 = ss12 = 0.0;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2) { mn1 += *X1; mn2 += *X2; }
        mn1 *= den; mn2 *= den; X1 -= L; X2 -= L;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            xx1 = *X1 - mn1; xx1 *= xx1;
            xx2 = *X2 - mn2; xx2 *= xx2;
            ss1 += xx1; ss2 += xx2; ss12 += xx1*xx2;
        }
        *Y = (double)L * ss12 / (ss1*ss2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? L : 0u, J2 = (L==N2) ? L : 0u;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                mn1 = mn2 = ss1 = ss2 = ss12 = 0.0;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2) { mn1 += *X1; mn2 += *X2; }
                mn1 *= den; mn2 *= den;
                X1 -= L; X2 -= L;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    xx1 = *X1 - mn1; xx1 *= xx1;
                    xx2 = *X2 - mn2; xx2 *= xx2;
                    ss1 += xx1; ss2 += xx2; ss12 += xx1*xx2;
                }
                *Y = (double)L * ss12 / (ss1*ss2);
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : B*(L-1u), I2 = (L==N2) ? 0u : B*(L-1u);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    mn1 = mn2 = ss1 = ss2 = ss12 = 0.0;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2) { mn1 += *X1; mn2 += *X2; }
                    mn1 *= den; mn2 *= den;
                    X1 -= L*K1; X2 -= L*K2;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        xx1 = *X1 - mn1; xx1 *= xx1;
                        xx2 = *X2 - mn2; xx2 *= xx2;
                        ss1 += xx1; ss2 += xx2; ss12 += xx1*xx2;
                    }
                    *Y = (double)L * ss12 / (ss1*ss2);
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
