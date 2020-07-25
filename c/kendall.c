//Vecs2scalar operation for 2 inputs X1 and X2.
//Kendall rank correlation coefficient for each pair of vectors.

#include <stdio.h>
#include <cblas.h>
#include <lapacke.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int kendall_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int kendall_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);


int kendall_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in kendall_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in kendall_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float den = 2.0f / (L*(L-1));
    int s1, s2, ssm = 0;

    if (N==0) {}
    else if (L==1)
    {
        const float o = 1.0f;
        cblas_scopy((int)N,&o,0,Y,1);
    }
    else if (L==N)
    {
        for (size_t l2=1; l2<L; l2++)
        {
            for (size_t l1=0; l1<l2; l1++)
            {
                s1 = (X1[l1]>X1[l2]) ? 1 : (X1[l1]<X1[l2]) ? -1 : 0;
                s2 = (X2[l1]>X2[l2]) ? 1 : (X2[l1]<X2[l2]) ? -1 : 0;
                ssm += s1 * s2;
            }
        }
        *Y = den * ssm;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? 0 : L, J2 = (L==N2) ? 0 : L;
            for (size_t v=0; v<V; v++, X1+=J1, X2+=J2)
            {
                ssm = 0;
                for (size_t l2=1; l2<L; l2++)
                {
                    for (size_t l1=0; l1<l2; l1++)
                    {
                        s1 = (X1[l1]>X1[l2]) ? 1 : (X1[l1]<X1[l2]) ? -1 : 0;
                        s2 = (X2[l1]>X2[l2]) ? 1 : (X2[l1]<X2[l2]) ? -1 : 0;
                        ssm += s1 * s2;
                    }
                }
                *Y++ = den * ssm;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            const size_t K1 = (L==N1) ? 1 : K, K2 = (L==N2) ? 1 : K;
            const size_t I1 = (L==N1) ? 0 : B*(L-1), I2 = (L==N2) ? 0 : B*(L-1);
            for (size_t g=0; g<G; g++, X1+=I1, X2+=I2)
            {
                for (size_t b=0; b<B; b++, X1+=J1, X2+=J2)
                {
                    ssm = 0;
                    for (size_t l2=1; l2<L; l2++)
                    {
                        for (size_t l1=0; l1<l2; l1++)
                        {
                            s1 = (X1[l1*K1]>X1[l2*K1]) ? 1 : (X1[l1*K1]<X1[l2*K1]) ? -1 : 0;
                            s2 = (X2[l1*K2]>X2[l2*K2]) ? 1 : (X2[l1*K2]<X2[l2*K2]) ? -1 : 0;
                            ssm += s1 * s2;
                        }
                    }
                    *Y++ = den * ssm;
                }
            }
        }
    }

    return 0;
}


int kendall_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in kendall_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in kendall_d: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double den = 2.0 / (L*(L-1));
    int s1, s2, ssm = 0;

    if (N==0) {}
    else if (L==1)
    {
        const double o = 1.0;
        cblas_dcopy((int)N,&o,0,Y,1);
    }
    else if (L==N)
    {
        for (size_t l2=1; l2<L; l2++)
        {
            for (size_t l1=0; l1<l2; l1++)
            {
                s1 = (X1[l1]>X1[l2]) ? 1 : (X1[l1]<X1[l2]) ? -1 : 0;
                s2 = (X2[l1]>X2[l2]) ? 1 : (X2[l1]<X2[l2]) ? -1 : 0;
                ssm += s1 * s2;
            }
        }
        *Y = den * ssm;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? 0 : L, J2 = (L==N2) ? 0 : L;
            for (size_t v=0; v<V; v++, X1+=J1, X2+=J2)
            {
                ssm = 0;
                for (size_t l2=1; l2<L; l2++)
                {
                    for (size_t l1=0; l1<l2; l1++)
                    {
                        s1 = (X1[l1]>X1[l2]) ? 1 : (X1[l1]<X1[l2]) ? -1 : 0;
                        s2 = (X2[l1]>X2[l2]) ? 1 : (X2[l1]<X2[l2]) ? -1 : 0;
                        ssm += s1 * s2;
                    }
                }
                *Y++ = den * ssm;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            const size_t K1 = (L==N1) ? 1 : K, K2 = (L==N2) ? 1 : K;
            const size_t I1 = (L==N1) ? 0 : B*(L-1), I2 = (L==N2) ? 0 : B*(L-1);
            for (size_t g=0; g<G; g++, X1+=I1, X2+=I2)
            {
                for (size_t b=0; b<B; b++, X1+=J1, X2+=J2)
                {
                    ssm = 0;
                    for (size_t l2=1; l2<L; l2++)
                    {
                        for (size_t l1=0; l1<l2; l1++)
                        {
                            s1 = (X1[l1*K1]>X1[l2*K1]) ? 1 : (X1[l1*K1]<X1[l2*K1]) ? -1 : 0;
                            s2 = (X2[l1*K2]>X2[l2*K2]) ? 1 : (X2[l1*K2]<X2[l2*K2]) ? -1 : 0;
                            ssm += s1 * s2;
                        }
                    }
                    *Y++ = den * ssm;
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
