//Vecs2scalar operation for 2 inputs X1 and X2.
//Correlation for each pair of vectors.
//This is the same as cov, except divides by std deviations of X1, X2.

//This is the Pearson product-moment correlation coefficient,
//a.k.a., the Pearson correlation coefficient, or Pearson's r,
//or the bivariate correlation, and is in [-1 1].

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

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
    if (dim>3) { fprintf(stderr,"error in corr_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in corr_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float ni = 1.0f/L;
    float x1, x2, mn1 = 0.0f, mn2 = 0.0f, sd1 = 0.0f, sd2 = 0.0f;

    if (N==0) {}
    else if (L==1)
    {
        const float o = 1.0f;
        cblas_scopy((int)N,&o,0,Y,1);
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; l++) { mn1 += *X1++; mn2 += *X2++; }
        mn1 *= ni; mn2 *= ni; X1 -= L; X2 -= L; *Y = 0.0f;
        for (size_t l=0; l<L; l++, X1++, X2++)
        {
            x1 = *X1 - mn1; x2 = *X2 - mn2;
            sd1 += x1*x1; sd2 += x2*x2; *Y += x1*x2;
        }
        *Y /= sqrtf(sd1*sd2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? L : 0, J2 = (L==N2) ? L : 0;
            for (size_t v=0; v<V; v++, X1-=J1, X2-=J2)
            {
                *Y = mn1 = mn2 = sd1 = sd2 = 0.0f;
                for (size_t l=0; l<L; l++) { mn1 += *X1++; mn2 += *X2++; }
                mn1 *= ni; mn2 *= ni; X1 -= L; X2 -= L;
                for (size_t l=0; l<L; l++, X1++, X2++)
                {
                    x1 = *X1 - mn1; x2 = *X2 - mn2;
                    sd1 += x1*x1; sd2 += x2*x2; *Y += x1*x2;
                }
                *Y++ /= sqrtf(sd1*sd2);
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            const size_t K1 = (L==N1) ? 1 : K, K2 = (L==N2) ? 1 : K;
            const size_t I1 = (L==N1) ? 0 : B*(L-1), I2 = (L==N2) ? 0 : B*(L-1);
            for (size_t g=0; g<G; g++, X1+=I1, X2+=I2)
            {
                for (size_t b=0; b<B; b++, X1-=L*K1-J1, X2-=L*K2-J2)
                {
                    *Y = mn1 = mn2 = sd1 = sd2 = 0.0f;
                    for (size_t l=0; l<L; l++, X1+=K1, X2+=K2) { mn1 += *X1; mn2 += *X2; }
                    mn1 *= ni; mn2 *= ni; X1 -= L*K1; X2 -= L*K2;
                    for (size_t l=0; l<L; l++, X1+=K1, X2+=K2)
                    {
                        x1 = *X1 - mn1; x2 = *X2 - mn2;
                        sd1 += x1*x1; sd2 += x2*x2; *Y += x1*x2;
                    }
                    *Y++ /= sqrtf(sd1*sd2);
                }
            }
        }
    }

    return 0;
}


int corr_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in corr_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in corr_d: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double ni = 1.0/L;
    double x1, x2, mn1 = 0.0, mn2 = 0.0, sd1 = 0.0, sd2 = 0.0;

    if (N==0) {}
    else if (L==1)
    {
        const double o = 1.0;
        cblas_dcopy((int)N,&o,0,Y,1);
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; l++) { mn1 += *X1++; mn2 += *X2++; }
        mn1 *= ni; mn2 *= ni; X1 -= L; X2 -= L; *Y = 0.0;
        for (size_t l=0; l<L; l++, X1++, X2++)
        {
            x1 = *X1 - mn1; x2 = *X2 - mn2;
            sd1 += x1*x1; sd2 += x2*x2; *Y += x1*x2;
        }
        *Y /= sqrt(sd1*sd2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? L : 0, J2 = (L==N2) ? L : 0;
            for (size_t v=0; v<V; v++, X1-=J1, X2-=J2)
            {
                *Y = mn1 = mn2 = sd1 = sd2 = 0.0;
                for (size_t l=0; l<L; l++) { mn1 += *X1++; mn2 += *X2++; }
                mn1 *= ni; mn2 *= ni; X1 -= L; X2 -= L;
                for (size_t l=0; l<L; l++, X1++, X2++)
                {
                    x1 = *X1 - mn1; x2 = *X2 - mn2;
                    sd1 += x1*x1; sd2 += x2*x2; *Y += x1*x2;
                }
                *Y++ /= sqrt(sd1*sd2);
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            const size_t K1 = (L==N1) ? 1 : K, K2 = (L==N2) ? 1 : K;
            const size_t I1 = (L==N1) ? 0 : B*(L-1), I2 = (L==N2) ? 0 : B*(L-1);
            for (size_t g=0; g<G; g++, X1+=I1, X2+=I2)
            {
                for (size_t b=0; b<B; b++, X1-=L*K1-J1, X2-=L*K2-J2)
                {
                    *Y = mn1 = mn2 = sd1 = sd2 = 0.0;
                    for (size_t l=0; l<L; l++, X1+=K1, X2+=K2) { mn1 += *X1; mn2 += *X2; }
                    mn1 *= ni; mn2 *= ni; X1 -= L*K1; X2 -= L*K2;
                    for (size_t l=0; l<L; l++, X1+=K1, X2+=K2)
                    {
                        x1 = *X1 - mn1; x2 = *X2 - mn2;
                        sd1 += x1*x1; sd2 += x2*x2; *Y += x1*x2;
                    }
                    *Y++ /= sqrt(sd1*sd2);
                }
            }
        }
    }

    return 0;
}


int corr_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in corr_c: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in corr_c: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float ni = 1.0f/L;
    float x1r, x1i, x2r, x2i, mn1r, mn1i, mn2r, mn2i, sd1, sd2, den;

    if (N==0) {}
    else if (L==1)
    {
        const float o = 1.0f;
        cblas_scopy((int)N,&o,0,Y,2);
    }
    else if (L==N)
    {
        mn1r = mn1i = mn2r = mn2i = sd1 = sd2 = 0.0f;
        for (size_t l=0; l<L; l++)
        {
            mn1r += *X1++; mn1i += *X1++;
            mn2r += *X2++; mn2i += *X2++;
        }
        mn1r *= ni; mn1i *= ni;
        mn2r *= ni; mn2i *= ni;
        X1 -= 2*L; X2 -= 2*L;
        *Y++ = 0.0f; *Y-- = 0.0f;
        for (size_t l=0; l<L; l++)
        {
            x1r = *X1++ - mn1r; x1i = *X1++ - mn1i;
            x2r = *X2++ - mn2r; x2i = *X2++ - mn2i;
            sd1 += x1r*x1r + x1i*x1i;
            sd2 += x2r*x2r + x2i*x2i;
            *Y++ += x1r*x2r + x1i*x2i;
            *Y-- -= x1r*x2i - x1i*x2r;
        }
        den = sqrtf(sd1*sd2);
        *Y++ /= den; *Y /= den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? 2*L : 0, J2 = (L==N2) ? 2*L : 0;
            for (size_t v=0; v<V; v++, X1-=J1, X2-=J2)
            {
                mn1r = mn1i = mn2r = mn2i = sd1 = sd2 = 0.0f;
                for (size_t l=0; l<L; l++)
                {
                    mn1r += *X1++; mn1i += *X1++;
                    mn2r += *X2++; mn2i += *X2++;
                }
                mn1r *= ni; mn1i *= ni;
                mn2r *= ni; mn2i *= ni;
                X1 -= 2*L; X2 -= 2*L;
                *Y++ = 0.0f; *Y-- = 0.0f;
                for (size_t l=0; l<L; l++)
                {
                    x1r = *X1++ - mn1r; x1i = *X1++ - mn1i;
                    x2r = *X2++ - mn2r; x2i = *X2++ - mn2i;
                    sd1 += x1r*x1r + x1i*x1i;
                    sd2 += x2r*x2r + x2i*x2i;
                    *Y++ += x1r*x2r + x1i*x2i;
                    *Y-- -= x1r*x2i - x1i*x2r;
                }
                den = sqrtf(sd1*sd2);
                *Y++ /= den; *Y++ /= den;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 2, J2 = (L==N2) ? 0 : 2;
            const size_t K1 = (L==N1) ? 2 : 2*K, K2 = (L==N2) ? 2 : 2*K;
            const size_t I1 = (L==N1) ? 0 : 2*B*(L-1), I2 = (L==N2) ? 0 : 2*B*(L-1);
            for (size_t g=0; g<G; g++, X1+=I1, X2+=I2)
            {
                for (size_t b=0; b<B; b++, X1-=L*K1-J1, X2-=L*K2-J2)
                {
                    mn1r = mn1i = mn2r = mn2i = sd1 = sd2 = 0.0f;
                    for (size_t l=0; l<L; l++, X1+=K1-1, X2+=K2-1)
                    {
                        mn1r += *X1++; mn1i += *X1;
                        mn2r += *X2++; mn2i += *X2;
                    }
                    mn1r *= ni; mn1i *= ni;
                    mn2r *= ni; mn2i *= ni;
                    X1 -= L*K1; X2 -= L*K2;
                    *Y++ = 0.0f; *Y-- = 0.0f;
                    for (size_t l=0; l<L; l++, X1+=K1-1, X2+=K2-1)
                    {
                        x1r = *X1++ - mn1r; x1i = *X1 - mn1i;
                        x2r = *X2++ - mn2r; x2i = *X2 - mn2i;
                        sd1 += x1r*x1r + x1i*x1i;
                        sd2 += x2r*x2r + x2i*x2i;
                        *Y++ += x1r*x2r + x1i*x2i;
                        *Y-- -= x1r*x2i - x1i*x2r;
                    }
                    den = sqrtf(sd1*sd2);
                    *Y++ /= den; *Y++ /= den;
                }
            }
        }
    }

    return 0;
}


int corr_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in corr_z: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in corr_z: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double ni = 1.0/L;
    double x1r, x1i, x2r, x2i, mn1r, mn1i, mn2r, mn2i, sd1, sd2, den;

    if (N==0) {}
    else if (L==1)
    {
        const double o = 1.0;
        cblas_dcopy((int)N,&o,0,Y,2);
    }
    else if (L==N)
    {
        mn1r = mn1i = mn2r = mn2i = sd1 = sd2 = 0.0;
        for (size_t l=0; l<L; l++)
        {
            mn1r += *X1++; mn1i += *X1++;
            mn2r += *X2++; mn2i += *X2++;
        }
        mn1r *= ni; mn1i *= ni;
        mn2r *= ni; mn2i *= ni;
        X1 -= 2*L; X2 -= 2*L;
        *Y++ = 0.0; *Y-- = 0.0;
        for (size_t l=0; l<L; l++)
        {
            x1r = *X1++ - mn1r; x1i = *X1++ - mn1i;
            x2r = *X2++ - mn2r; x2i = *X2++ - mn2i;
            sd1 += x1r*x1r + x1i*x1i;
            sd2 += x2r*x2r + x2i*x2i;
            *Y++ += x1r*x2r + x1i*x2i;
            *Y-- -= x1r*x2i - x1i*x2r;
        }
        den = sqrt(sd1*sd2);
        *Y++ /= den; *Y /= den;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? 2*L : 0, J2 = (L==N2) ? 2*L : 0;
            for (size_t v=0; v<V; v++, X1-=J1, X2-=J2)
            {
                mn1r = mn1i = mn2r = mn2i = sd1 = sd2 = 0.0;
                for (size_t l=0; l<L; l++)
                {
                    mn1r += *X1++; mn1i += *X1++;
                    mn2r += *X2++; mn2i += *X2++;
                }
                mn1r *= ni; mn1i *= ni;
                mn2r *= ni; mn2i *= ni;
                X1 -= 2*L; X2 -= 2*L;
                *Y++ = 0.0; *Y-- = 0.0;
                for (size_t l=0; l<L; l++)
                {
                    x1r = *X1++ - mn1r; x1i = *X1++ - mn1i;
                    x2r = *X2++ - mn2r; x2i = *X2++ - mn2i;
                    sd1 += x1r*x1r + x1i*x1i;
                    sd2 += x2r*x2r + x2i*x2i;
                    *Y++ += x1r*x2r + x1i*x2i;
                    *Y-- -= x1r*x2i - x1i*x2r;
                }
                den = sqrt(sd1*sd2);
                *Y++ /= den; *Y++ /= den;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 2, J2 = (L==N2) ? 0 : 2;
            const size_t K1 = (L==N1) ? 2 : 2*K, K2 = (L==N2) ? 2 : 2*K;
            const size_t I1 = (L==N1) ? 0 : 2*B*(L-1), I2 = (L==N2) ? 0 : 2*B*(L-1);
            for (size_t g=0; g<G; g++, X1+=I1, X2+=I2)
            {
                for (size_t b=0; b<B; b++, X1-=L*K1-J1, X2-=L*K2-J2)
                {
                    mn1r = mn1i = mn2r = mn2i = sd1 = sd2 = 0.0;
                    for (size_t l=0; l<L; l++, X1+=K1-1, X2+=K2-1)
                    {
                        mn1r += *X1++; mn1i += *X1;
                        mn2r += *X2++; mn2i += *X2;
                    }
                    mn1r *= ni; mn1i *= ni;
                    mn2r *= ni; mn2i *= ni;
                    X1 -= L*K1; X2 -= L*K2;
                    *Y++ = 0.0; *Y-- = 0.0;
                    for (size_t l=0; l<L; l++, X1+=K1-1, X2+=K2-1)
                    {
                        x1r = *X1++ - mn1r; x1i = *X1 - mn1i;
                        x2r = *X2++ - mn2r; x2i = *X2 - mn2i;
                        sd1 += x1r*x1r + x1i*x1i;
                        sd2 += x2r*x2r + x2i*x2i;
                        *Y++ += x1r*x2r + x1i*x2i;
                        *Y-- -= x1r*x2i - x1i*x2r;
                    }
                    den = sqrt(sd1*sd2);
                    *Y++ /= den; *Y++ /= den;
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
