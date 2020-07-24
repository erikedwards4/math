//Vecs2scalar operation for 2 inputs X1 and X2.
//Covariance for each pair of vectors.
//This is the normalized dot product after removing the means.
//The normalization is either 1/L or 1/(L-1), depending on biased.

//For complex inputs, this uses the conjugated dot product.

//Note that this is a function for vec operations only.
//This does NOT compute a covariance matrix for multi-channel input.
//This computes pairwise covariance taken as a vector similarity measure.

#include <stdio.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cov_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const char biased);
int cov_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const char biased);
int cov_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const char biased);
int cov_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const char biased);


int cov_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in cov_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in cov_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float ni = 1.0f/L;
    const float den = (biased) ? ni : 1.0f/(L-1);
    float mn1 = 0.0f, mn2 = 0.0f;

    if (N==0) {}
    else if (L==N)
    {
        for (size_t l=0; l<L; l++) { mn1 += *X1++; mn2 += *X2++; }
        mn1 *= ni; mn2 *= ni; X1 -= L; X2 -= L; *Y = 0.0f;
        for (size_t l=0; l<L; l++, X1++, X2++) { *Y += (*X1-mn1) * (*X2-mn2); }
        *Y *= den;
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
                mn1 = mn2 = *Y = 0.0f;
                for (size_t l=0; l<L; l++) { mn1 += *X1++; mn2 += *X2++; }
                mn1 *= ni; mn2 *= ni; X1 -= L; X2 -= L;
                for (size_t l=0; l<L; l++, X1++, X2++) { *Y += (*X1-mn1) * (*X2-mn2); }
                *Y++ *= den;
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
                    mn1 = mn2 = *Y = 0.0f;
                    for (size_t l=0; l<L; l++, X1+=K1, X2+=K2) { mn1 += *X1; mn2 += *X2; }
                    mn1 *= ni; mn2 *= ni; X1 -= L*K1; X2 -= L*K2;
                    for (size_t l=0; l<L; l++, X1+=K1, X2+=K2) { *Y += (*X1-mn1) * (*X2-mn2); }
                    *Y++ *= den;
                }
            }
        }
    }

    return 0;
}


int cov_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in cov_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in cov_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double ni = 1.0/L;
    const double den = (biased) ? ni : 1.0/(L-1);
    double mn1 = 0.0, mn2 = 0.0;

    if (N==0) {}
    else if (L==N)
    {
        for (size_t l=0; l<L; l++) { mn1 += *X1++; mn2 += *X2++; }
        mn1 *= ni; mn2 *= ni; X1 -= L; X2 -= L; *Y = 0.0;
        for (size_t l=0; l<L; l++, X1++, X2++) { *Y += (*X1-mn1) * (*X2-mn2); }
        *Y *= den;
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
                mn1 = mn2 = *Y = 0.0;
                for (size_t l=0; l<L; l++) { mn1 += *X1++; mn2 += *X2++; }
                mn1 *= ni; mn2 *= ni; X1 -= L; X2 -= L;
                for (size_t l=0; l<L; l++, X1++, X2++) { *Y += (*X1-mn1) * (*X2-mn2); }
                *Y++ *= den;
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
                    mn1 = mn2 = *Y = 0.0;
                    for (size_t l=0; l<L; l++, X1+=K1, X2+=K2) { mn1 += *X1; mn2 += *X2; }
                    mn1 *= ni; mn2 *= ni; X1 -= L*K1; X2 -= L*K2;
                    for (size_t l=0; l<L; l++, X1+=K1, X2+=K2) { *Y += (*X1-mn1) * (*X2-mn2); }
                    *Y++ *= den;
                }
            }
        }
    }

    return 0;
}


int cov_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in cov_c: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in cov_c: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float ni = 1.0f/L;
    const float den = (biased) ? ni : 1.0f/(L-1);
    float mn1r, mn1i, mn2r, mn2i, x1r, x1i, x2r, x2i;

    if (N==0) {}
    else if (L==N)
    {
        mn1r = mn1i = mn2r = mn2i = 0.0f;
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
            *Y++ += x1r*x2r + x1i*x2i;
            *Y-- -= x1r*x2i - x1i*x2r;
        }
        *Y++ *= den; *Y *= den;
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
                mn1r = mn1i = mn2r = mn2i = 0.0f;
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
                    *Y++ += x1r*x2r + x1i*x2i;
                    *Y-- -= x1r*x2i - x1i*x2r;
                }
                *Y++ *= den; *Y++ *= den;
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
                    mn1r = mn1i = mn2r = mn2i = 0.0f;
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
                        *Y++ += x1r*x2r + x1i*x2i;
                        *Y-- -= x1r*x2i - x1i*x2r;
                    }
                    *Y++ *= den; *Y++ *= den;
                }
            }
        }
    }

    return 0;
}


int cov_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const char biased)
{
    if (dim>3) { fprintf(stderr,"error in cov_z: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in cov_z: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double ni = 1.0/L;
    const double den = (biased) ? ni : 1.0/(L-1);
    double mn1r, mn1i, mn2r, mn2i, x1r, x1i, x2r, x2i;

    if (N==0) {}
    else if (L==N)
    {
        mn1r = mn1i = mn2r = mn2i = 0.0;
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
            *Y++ += x1r*x2r + x1i*x2i;
            *Y-- -= x1r*x2i - x1i*x2r;
        }
        *Y++ *= den; *Y *= den;
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
                mn1r = mn1i = mn2r = mn2i = 0.0;
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
                    *Y++ += x1r*x2r + x1i*x2i;
                    *Y-- -= x1r*x2i - x1i*x2r;
                }
                *Y++ *= den; *Y++ *= den;
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
                    mn1r = mn1i = mn2r = mn2i = 0.0;
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
                        *Y++ += x1r*x2r + x1i*x2i;
                        *Y-- -= x1r*x2i - x1i*x2r;
                    }
                    *Y++ *= den; *Y++ *= den;
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
