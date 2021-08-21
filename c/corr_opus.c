//Vecs2scalar operation for 2 inputs X1 and X2.
//Correlation for each pair of vectors,
//but using C_opus as defined by Eq. (2) of Vos et al. [2013].
//This is the same as cov, except divides by
//the mean of sum(X1.^2) and sum(X2.^2), i.e.:

// C_opus = 2.0 * x1'*x2 / (x1'*x1 + x2'*x2)

//"This correlation measures similarity not just in shape, but also in scale.
// Two vectors with very different energies will have a lower correlation,
// similar to frequency-domain pitch estimators."

//Vos K, Sørensen KV, Jensen SS, Valin JM. 2013. Voice coding with Opus. Proc AES Convention. Audio Engineering Society: paper 8941.

//Unfortunately, they do not specify if they subtract the means from X1 and X2 first.
//Since it is a correlation, and they show the usual corr also, I assume that they do.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int corr_opus_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
int corr_opus_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
int corr_opus_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);
int corr_opus_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim);


int corr_opus_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in corr_opus_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in corr_opus_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float den = 1.0f/(float)L;
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
        *Y = 2.0f * sm2 / (sd1+sd2);
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
                mn1 = mn2 = sd1 = sd2 = sm2 = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2) { mn1 += *X1; mn2 += *X2; }
                mn1 *= den; mn2 *= den;
                X1 -= L; X2 -= L;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    x1 = *X1 - mn1; x2 = *X2 - mn2;
                    sd1 += x1*x1; sd2 += x2*x2; sm2 += x1*x2;
                }
                *Y = 2.0f * sm2 / (sd1+sd2);
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
                    mn1 = mn2 = sd1 = sd2 = sm2 = 0.0f;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2) { mn1 += *X1; mn2 += *X2; }
                    mn1 *= den; mn2 *= den;
                    X1 -= L*K1; X2 -= L*K2;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        x1 = *X1 - mn1; x2 = *X2 - mn2;
                        sd1 += x1*x1; sd2 += x2*x2; sm2 += x1*x2;
                    }
                    *Y = 2.0f * sm2 / (sd1+sd2);
                }
            }
        }
    }

    return 0;
}


int corr_opus_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in corr_opus_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in corr_opus_d: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double den = 1.0/(double)L;
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
        *Y = 2.0 * sm2 / (sd1+sd2);
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
                mn1 = mn2 = sd1 = sd2 = sm2 = 0.0;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2) { mn1 += *X1; mn2 += *X2; }
                mn1 *= den; mn2 *= den;
                X1 -= L; X2 -= L;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    x1 = *X1 - mn1; x2 = *X2 - mn2;
                    sd1 += x1*x1; sd2 += x2*x2; sm2 += x1*x2;
                }
                *Y = 2.0 * sm2 / (sd1+sd2);
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
                    mn1 = mn2 = sd1 = sd2 = sm2 = 0.0;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2) { mn1 += *X1; mn2 += *X2; }
                    mn1 *= den; mn2 *= den;
                    X1 -= L*K1; X2 -= L*K2;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        x1 = *X1 - mn1; x2 = *X2 - mn2;
                        sd1 += x1*x1; sd2 += x2*x2; sm2 += x1*x2;
                    }
                    *Y = 2.0 * sm2 / (sd1+sd2);
                }
            }
        }
    }

    return 0;
}


int corr_opus_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in corr_opus_c: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in corr_opus_c: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float den = 1.0f/(float)L;
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
        den2 = 0.5f * (sd1+sd2);
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
            for (size_t v=V; v>0u; --v, X1-=J1, X2-=J2)
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
                den2 = 0.5f * (sd1+sd2);
                *Y = yr / den2; *++Y = yi / den2;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 2u, J2 = (L==N2) ? 0u : 2u;
            const size_t K1 = (L==N1) ? 2u : 2u*K, K2 = (L==N2) ? 2u : 2u*K;
            const size_t I1 = (L==N1) ? 0u : 2u*B*(L-1u), I2 = (L==N2) ? 0u : 2u*B*(L-1u);
            for (size_t g=G; g>0u; --g, X1+=I1, X2+=I2)
            {
                for (size_t b=B; b>0u; --b, X1-=L*K1-J1, X2-=L*K2-J2)
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
                    den2 = 0.5f * (sd1+sd2);
                    *Y = yr / den2; *++Y = yi / den2;
                }
            }
        }
    }

    return 0;
}


int corr_opus_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in corr_opus_z: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in corr_opus_z: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double den = 1.0/(double)L;
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
        den2 = 0.5 * (sd1+sd2);
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
            for (size_t v=V; v>0u; --v, X1-=J1, X2-=J2)
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
                den2 = 0.5 * (sd1+sd2);
                *Y = yr / den2; *++Y = yi / den2;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 2u, J2 = (L==N2) ? 0u : 2u;
            const size_t K1 = (L==N1) ? 2u : 2u*K, K2 = (L==N2) ? 2u : 2u*K;
            const size_t I1 = (L==N1) ? 0u : 2u*B*(L-1u), I2 = (L==N2) ? 0u : 2u*B*(L-1u);
            for (size_t g=G; g>0u; --g, X1+=I1, X2+=I2)
            {
                for (size_t b=B; b>0u; --b, X1-=L*K1-J1, X2-=L*K2-J2)
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
                    den2 = 0.5 * (sd1+sd2);
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
