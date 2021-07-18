//Vecs2scalar operation for 2 inputs X1 and X2.
//Correlation for each pair of vectors.
//This is the same as cov, except divides by std deviations of X1, X2.

//This is the Pearson product-moment correlation coefficient,
//a.k.a., the Pearson correlation coefficient, or Pearson's r,
//or the bivariate correlation, and is in [-1 1].

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int corr_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int corr_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int corr_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int corr_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);


int corr_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in corr_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in corr_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float den = 1.0f/L;
    float x1, x2, mn1 = 0.0f, mn2 = 0.0f, sd1 = 0.0f, sd2 = 0.0f, sm2 = 0.0f;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 1.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X1, ++X2) { mn1 += *X1; mn2 += *X2; }
        mn1 *= den; mn2 *= den;
        X1 -= L; X2 -= L;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            x1 = *X1 - mn1; x2 = *X2 - mn2;
            sd1 += x1*x1; sd2 += x2*x2; sm2 += x1*x2;
        }
        *Y = sm2 / sqrtf(sd1*sd2);
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
                mn1 = mn2 = sd1 = sd2 = sm2 = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2) { mn1 += *X1; mn2 += *X2; }
                mn1 *= den; mn2 *= den;
                X1 -= L; X2 -= L;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    x1 = *X1 - mn1; x2 = *X2 - mn2;
                    sd1 += x1*x1; sd2 += x2*x2; sm2 += x1*x2;
                }
                *Y = sm2 * sqrtf(sd1*sd2);
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
                    mn1 = mn2 = sd1 = sd2 = sm2 = 0.0f;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2) { mn1 += *X1; mn2 += *X2; }
                    mn1 *= den; mn2 *= den;
                    X1 -= L*K1; X2 -= L*K2;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        x1 = *X1 - mn1; x2 = *X2 - mn2;
                        sd1 += x1*x1; sd2 += x2*x2; sm2 += x1*x2;
                    }
                    *Y = sm2 / sqrtf(sd1*sd2);
                }
            }
        }
    }

    return 0;
}


int corr_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in corr_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in corr_d: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double den = 1.0/L;
    double x1, x2, mn1 = 0.0, mn2 = 0.0, sd1 = 0.0, sd2 = 0.0, sm2 = 0.0;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 1.0; }
    }
    else if (L==N)
    {
        for (size_t l=0u; l<L; ++l, ++X1, ++X2) { mn1 += *X1; mn2 += *X2; }
        mn1 *= den; mn2 *= den;
        X1 -= L; X2 -= L;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            x1 = *X1 - mn1; x2 = *X2 - mn2;
            sd1 += x1*x1; sd2 += x2*x2; sm2 += x1*x2;
        }
        *Y = sm2 / sqrt(sd1*sd2);
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
                mn1 = mn2 = sd1 = sd2 = sm2 = 0.0;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2) { mn1 += *X1; mn2 += *X2; }
                mn1 *= den; mn2 *= den;
                X1 -= L; X2 -= L;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    x1 = *X1 - mn1; x2 = *X2 - mn2;
                    sd1 += x1*x1; sd2 += x2*x2; sm2 += x1*x2;
                }
                *Y = sm2 * sqrt(sd1*sd2);
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
                    mn1 = mn2 = sd1 = sd2 = sm2 = 0.0;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2) { mn1 += *X1; mn2 += *X2; }
                    mn1 *= den; mn2 *= den;
                    X1 -= L*K1; X2 -= L*K2;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        x1 = *X1 - mn1; x2 = *X2 - mn2;
                        sd1 += x1*x1; sd2 += x2*x2; sm2 += x1*x2;
                    }
                    *Y = sm2 / sqrt(sd1*sd2);
                }
            }
        }
    }

    return 0;
}


int corr_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in corr_c: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in corr_c: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float den = 1.0f/L;
    float x1r, x1i, x2r, x2i, mn1r, mn1i, mn2r, mn2i, sd1, sd2, yr, yi, den2;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 1.0f; }
    }
    else if (L==N)
    {
        mn1r = mn1i = mn2r = mn2i = sd1 = sd2 = yr = yi = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            mn1r += *X1; mn1i += *++X1;
            mn2r += *X2; mn2i += *++X2;
        }
        mn1r *= den; mn1i *= den;
        mn2r *= den; mn2i *= den;
        X1 -= 2u*L; X2 -= 2u*L;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            x1r = *X1 - mn1r; x1i = *++X1 - mn1i;
            x2r = *X2 - mn2r; x2i = *++X2 - mn2i;
            sd1 += x1r*x1r + x1i*x1i;
            sd2 += x2r*x2r + x2i*x2i;
            yr += x1r*x2r + x1i*x2i;
            yi -= x1r*x2i - x1i*x2r;
        }
        den2 = sqrtf(sd1*sd2);
        *Y = yr / den2; *++Y = yi / den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? 2u*L : 0u, J2 = (L==N2) ? 2u*L : 0u;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2)
            {
                mn1r = mn1i = mn2r = mn2i = sd1 = sd2 = yr = yi = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    mn1r += *X1; mn1i += *++X1;
                    mn2r += *X2; mn2i += *++X2;
                }
                mn1r *= den; mn1i *= den;
                mn2r *= den; mn2i *= den;
                X1 -= 2u*L; X2 -= 2u*L;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2, ++Y)
                {
                    x1r = *X1 - mn1r; x1i = *++X1 - mn1i;
                    x2r = *X2 - mn2r; x2i = *++X2 - mn2i;
                    sd1 += x1r*x1r + x1i*x1i;
                    sd2 += x2r*x2r + x2i*x2i;
                    yr += x1r*x2r + x1i*x2i;
                    yi -= x1r*x2i - x1i*x2r;
                }
                den2 = sqrtf(sd1*sd2);
                *Y = yr / den2; *++Y = yi / den2;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 2u, J2 = (L==N2) ? 0u : 2u;
            const size_t K1 = (L==N1) ? 2u : 2u*K, K2 = (L==N2) ? 2u : 2u*K;
            const size_t I1 = (L==N1) ? 0u : 2u*B*(L-1u), I2 = (L==N2) ? 0u : 2u*B*(L-1u);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2)
                {
                    mn1r = mn1i = mn2r = mn2i = sd1 = sd2 = yr = yi = 0.0f;
                    for (size_t l=0u; l<L; ++l, X1+=K1-1u, X2+=K2-1u)
                    {
                        mn1r += *X1; mn1i += *++X1;
                        mn2r += *X2; mn2i += *++X2;
                    }
                    mn1r *= den; mn1i *= den;
                    mn2r *= den; mn2i *= den;
                    X1 -= L*K1; X2 -= L*K2;
                    for (size_t l=0u; l<L; ++l, X1+=K1-1u, X2+=K2-1u, ++Y)
                    {
                        x1r = *X1 - mn1r; x1i = *++X1 - mn1i;
                        x2r = *X2 - mn2r; x2i = *++X2 - mn2i;
                        sd1 += x1r*x1r + x1i*x1i;
                        sd2 += x2r*x2r + x2i*x2i;
                        yr += x1r*x2r + x1i*x2i;
                        yi -= x1r*x2i - x1i*x2r;
                    }
                    den2 = sqrtf(sd1*sd2);
                    *Y = yr / den2; *++Y = yi / den2;
                }
            }
        }
    }

    return 0;
}


int corr_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in corr_z: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in corr_z: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double den = 1.0/L;
    double x1r, x1i, x2r, x2i, mn1r, mn1i, mn2r, mn2i, sd1, sd2, yr, yi, den2;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 1.0; }
    }
    else if (L==N)
    {
        mn1r = mn1i = mn2r = mn2i = sd1 = sd2 = yr = yi = 0.0;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            mn1r += *X1; mn1i += *++X1;
            mn2r += *X2; mn2i += *++X2;
        }
        mn1r *= den; mn1i *= den;
        mn2r *= den; mn2i *= den;
        X1 -= 2u*L; X2 -= 2u*L;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            x1r = *X1 - mn1r; x1i = *++X1 - mn1i;
            x2r = *X2 - mn2r; x2i = *++X2 - mn2i;
            sd1 += x1r*x1r + x1i*x1i;
            sd2 += x2r*x2r + x2i*x2i;
            yr += x1r*x2r + x1i*x2i;
            yi -= x1r*x2i - x1i*x2r;
        }
        den2 = sqrt(sd1*sd2);
        *Y = yr / den2; *++Y = yi / den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? 2u*L : 0u, J2 = (L==N2) ? 2u*L : 0u;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2)
            {
                mn1r = mn1i = mn2r = mn2i = sd1 = sd2 = yr = yi = 0.0;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    mn1r += *X1; mn1i += *++X1;
                    mn2r += *X2; mn2i += *++X2;
                }
                mn1r *= den; mn1i *= den;
                mn2r *= den; mn2i *= den;
                X1 -= 2u*L; X2 -= 2u*L;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2, ++Y)
                {
                    x1r = *X1 - mn1r; x1i = *++X1 - mn1i;
                    x2r = *X2 - mn2r; x2i = *++X2 - mn2i;
                    sd1 += x1r*x1r + x1i*x1i;
                    sd2 += x2r*x2r + x2i*x2i;
                    yr += x1r*x2r + x1i*x2i;
                    yi -= x1r*x2i - x1i*x2r;
                }
                den2 = sqrt(sd1*sd2);
                *Y = yr / den2; *++Y = yi / den2;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 2u, J2 = (L==N2) ? 0u : 2u;
            const size_t K1 = (L==N1) ? 2u : 2u*K, K2 = (L==N2) ? 2u : 2u*K;
            const size_t I1 = (L==N1) ? 0u : 2u*B*(L-1u), I2 = (L==N2) ? 0u : 2u*B*(L-1u);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2)
                {
                    mn1r = mn1i = mn2r = mn2i = sd1 = sd2 = yr = yi = 0.0;
                    for (size_t l=0u; l<L; ++l, X1+=K1-1u, X2+=K2-1u)
                    {
                        mn1r += *X1; mn1i += *++X1;
                        mn2r += *X2; mn2i += *++X2;
                    }
                    mn1r *= den; mn1i *= den;
                    mn2r *= den; mn2i *= den;
                    X1 -= L*K1; X2 -= L*K2;
                    for (size_t l=0u; l<L; ++l, X1+=K1-1u, X2+=K2-1u, ++Y)
                    {
                        x1r = *X1 - mn1r; x1i = *++X1 - mn1i;
                        x2r = *X2 - mn2r; x2i = *++X2 - mn2i;
                        sd1 += x1r*x1r + x1i*x1i;
                        sd2 += x2r*x2r + x2i*x2i;
                        yr += x1r*x2r + x1i*x2i;
                        yi -= x1r*x2i - x1i*x2r;
                    }
                    den2 = sqrt(sd1*sd2);
                    *Y = yr / den2; *++Y = yi / den2;
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
