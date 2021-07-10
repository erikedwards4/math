//Vecs2scalar operation for 2 inputs X1 and X2.
//Cosine similarity for each pair of vectors.
//This is the same as corr, except no mean subtraction.
//This is the dot product with L2-norm normalization, and is in [-1 1].
//This is equal to the cosine of the angle between the vectors.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cos2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int cos2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int cos2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int cos2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);


int cos2_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cos2_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in cos2_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    float sm2, sd1, sd2;

    if (N==0u) {}
    else if (L==N)
    {
        sm2 = sd1 = sd2 = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            sd1 += *X1**X1; sd2 += *X2**X2;
            sm2 += *X1**X2;
        }
        *Y = sm2 / sqrtf(sd1*sd2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? L : 0, J2 = (L==N2) ? L : 0;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                sm2 = sd1 = sd2 = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    sd1 += *X1**X1; sd2 += *X2**X2;
                    sm2 += *X1**X2;
                }
                *Y = sm2 / sqrtf(sd1*sd2);
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            const size_t K1 = (L==N1) ? 1 : K, K2 = (L==N2) ? 1 : K;
            const size_t I1 = (L==N1) ? 0 : B*(L-1), I2 = (L==N2) ? 0 : B*(L-1);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    sm2 = sd1 = sd2 = 0.0f;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        sd1 += *X1**X1; sd2 += *X2**X2;
                        sm2 += *X1**X2;
                    }
                    *Y = sm2 / sqrtf(sd1*sd2);
                }
            }
        }
    }

    return 0;
}


int cos2_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cos2_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in cos2_d: vectors in X1 and X2 must have the same length\n"); return 1; }
    double sm2, sd1, sd2;

    if (N==0u) {}
    else if (L==N)
    {
        sm2 = sd1 = sd2 = 0.0;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            sd1 += *X1**X1; sd2 += *X2**X2;
            sm2 += *X1**X2;
        }
        *Y = sm2 / sqrt(sd1*sd2);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? L : 0, J2 = (L==N2) ? L : 0;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                sm2 = sd1 = sd2 = 0.0;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    sd1 += *X1**X1; sd2 += *X2**X2;
                    sm2 += *X1**X2;
                }
                *Y = sm2 / sqrt(sd1*sd2);
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            const size_t K1 = (L==N1) ? 1 : K, K2 = (L==N2) ? 1 : K;
            const size_t I1 = (L==N1) ? 0 : B*(L-1), I2 = (L==N2) ? 0 : B*(L-1);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    sm2 = sd1 = sd2 = 0.0;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        sd1 += *X1**X1; sd2 += *X2**X2;
                        sm2 += *X1**X2;
                    }
                    *Y = sm2 / sqrt(sd1*sd2);
                }
            }
        }
    }

    return 0;
}


int cos2_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cos2_c: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in cos2_c: vectors in X1 and X2 must have the same length\n"); return 1; }
    float x1r, x2r, sd1, sd2, yr, yi, den2;

    if (N==0u) {}
    else if (L==N)
    {
        sd1 = sd2 = yr = yi = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            x1r = *X1++; x2r = *X2++;
            sd1 += x1r*x1r + *X1**X1;
            sd2 += x2r*x2r + *X2**X2;
            yr += x1r*x2r + *X1**X2;
            yi -= x1r**X2 - *X1*x2r;
        }
        den2 = sqrtf(sd1*sd2);
        *Y = yr / den2; *++Y = yi / den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? 2*L : 0, J2 = (L==N2) ? 2*L : 0;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                sd1 = sd2 = yr = yi = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2, ++Y)
                {
                    x1r = *X1++; x2r = *X2++;
                    sd1 += x1r*x1r + *X1**X1;
                    sd2 += x2r*x2r + *X2**X2;
                    yr += x1r*x2r + *X1**X2;
                    yi -= x1r**X2 - *X1*x2r;
                }
                den2 = sqrtf(sd1*sd2);
                *Y = yr / den2; *++Y = yi / den2;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 2, J2 = (L==N2) ? 0 : 2;
            const size_t K1 = (L==N1) ? 2 : 2*K, K2 = (L==N2) ? 2 : 2*K;
            const size_t I1 = (L==N1) ? 0 : 2*B*(L-1), I2 = (L==N2) ? 0 : 2*B*(L-1);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2)
                {
                    sd1 = sd2 = yr = yi = 0.0f;
                    for (size_t l=0u; l<L; ++l, X1+=K1-1, X2+=K2-1, ++Y)
                    {
                        x1r = *X1++; x2r = *X2++;
                        sd1 += x1r*x1r + *X1**X1;
                        sd2 += x2r*x2r + *X2**X2;
                        yr += x1r*x2r + *X1**X2;
                        yi -= x1r**X2 - *X1*x2r;
                    }
                    den2 = sqrtf(sd1*sd2);
                    *Y = yr / den2; *++Y = yi / den2;
                }
            }
        }
    }

    return 0;
}


int cos2_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in cos2_z: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in cos2_z: vectors in X1 and X2 must have the same length\n"); return 1; }
    double x1r, x2r, sd1, sd2, yr, yi, den2;

    if (N==0u) {}
    else if (L==N)
    {
        sd1 = sd2 = yr = yi = 0.0;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            x1r = *X1++; x2r = *X2++;
            sd1 += x1r*x1r + *X1**X1;
            sd2 += x2r*x2r + *X2**X2;
            yr += x1r*x2r + *X1**X2;
            yi -= x1r**X2 - *X1*x2r;
        }
        den2 = sqrt(sd1*sd2);
        *Y = yr / den2; *++Y = yi / den2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? 2*L : 0, J2 = (L==N2) ? 2*L : 0;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                sd1 = sd2 = yr = yi = 0.0;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2, ++Y)
                {
                    x1r = *X1++; x2r = *X2++;
                    sd1 += x1r*x1r + *X1**X1;
                    sd2 += x2r*x2r + *X2**X2;
                    yr += x1r*x2r + *X1**X2;
                    yi -= x1r**X2 - *X1*x2r;
                }
                den2 = sqrt(sd1*sd2);
                *Y = yr / den2; *++Y = yi / den2;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 2, J2 = (L==N2) ? 0 : 2;
            const size_t K1 = (L==N1) ? 2 : 2*K, K2 = (L==N2) ? 2 : 2*K;
            const size_t I1 = (L==N1) ? 0 : 2*B*(L-1), I2 = (L==N2) ? 0 : 2*B*(L-1);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2)
                {
                    sd1 = sd2 = yr = yi = 0.0;
                    for (size_t l=0u; l<L; ++l, X1+=K1-1, X2+=K2-1, ++Y)
                    {
                        x1r = *X1++; x2r = *X2++;
                        sd1 += x1r*x1r + *X1**X1;
                        sd2 += x2r*x2r + *X2**X2;
                        yr += x1r*x2r + *X1**X2;
                        yi -= x1r**X2 - *X1*x2r;
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
