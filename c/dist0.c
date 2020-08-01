//Vecs2scalar operation for 2 inputs X1 and X2.
//Gets L0 distance between pairs of vecs in X1, X2.

//This is the 'Hamming distance',
//which is just the count of non-zero values in |X1-X2|,
//which is just the count of X1!=X2.
//However, for numerical purposes, I use |X1-X2|>thresh.

#include <stdio.h>
#include <float.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int dist0_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int dist0_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int dist0_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int dist0_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);


int dist0_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in dist0_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist0_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    float thresh = 2.0f * FLT_EPSILON;

    if (N==0) {}
    else if (L==N)
    {
        size_t cnt = 0;
        for (size_t l=0; l<L; ++l, ++X1, ++X2)
        {
            cnt += (fabsf(*X1-*X2)>thresh);
        }
        *Y = cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? L : 0, J2 = (L==N2) ? L : 0;
            for (size_t v=0; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                size_t cnt = 0;
                for (size_t l=0; l<L; ++l, ++X1, ++X2)
                {
                    //cnt += (*X1!=*X2);
                    cnt += (fabsf(*X1-*X2)>thresh);
                }
                *Y = cnt;
            }
        }
        else if (G==1)
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            for (size_t v=0; v<V; ++v, X1+=J1, X2+=J2, ++Y)
            {
                //*Y = (float)(*X1!=*X2);
                *Y = (float)(fabsf(*X1-*X2)>thresh);
            }
            X1 += 1-J1; X2 += 1-J2; Y -= V;
            for (size_t l=1; l<L; ++l, X1+=1-J1, X2+=1-J2, Y-=V)
            {
                for (size_t v=0; v<V; ++v, X1+=J1, X2+=J2, ++Y)
                {
                    //*Y += (float)(*X1!=*X2);
                    *Y += (float)(fabsf(*X1-*X2)>thresh);
                }
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            const size_t K1 = (L==N1) ? 1 : K, K2 = (L==N2) ? 1 : K;
            const size_t I1 = (L==N1) ? 0 : B*(L-1), I2 = (L==N2) ? 0 : B*(L-1);
            for (size_t g=0; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2)
                {
                    size_t cnt = 0;
                    for (size_t l=0; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        //cnt += (*X1!=*X2);
                        cnt += (fabsf(*X1-*X2)>thresh);
                    }
                    *Y = cnt;
                }
            }
        }
    }

    return 0;
}


int dist0_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in dist0_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist0_d: vectors in X1 and X2 must have the same length\n"); return 1; }
    double thresh = 2.0 * DBL_EPSILON;

    if (N==0) {}
    else if (L==N)
    {
        size_t cnt = 0;
        for (size_t l=0; l<L; ++l, ++X1, ++X2)
        {
            cnt += (fabs(*X1-*X2)>thresh);
        }
        *Y = cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? L : 0, J2 = (L==N2) ? L : 0;
            for (size_t v=0; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                size_t cnt = 0;
                for (size_t l=0; l<L; ++l, ++X1, ++X2)
                {
                    cnt += (fabs(*X1-*X2)>thresh);
                }
                *Y = cnt;
            }
        }
        else if (G==1)
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            for (size_t v=0; v<V; ++v, X1+=J1, X2+=J2, ++Y)
            {
                *Y = (double)(fabs(*X1-*X2)>thresh);
            }
            X1 += 1-J1; X2 += 1-J2; Y -= V;
            for (size_t l=1; l<L; ++l, X1+=1-J1, X2+=1-J2, Y-=V)
            {
                for (size_t v=0; v<V; ++v, X1+=J1, X2+=J2, ++Y)
                {
                    *Y += (double)(fabs(*X1-*X2)>thresh);
                }
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            const size_t K1 = (L==N1) ? 1 : K, K2 = (L==N2) ? 1 : K;
            const size_t I1 = (L==N1) ? 0 : B*(L-1), I2 = (L==N2) ? 0 : B*(L-1);
            for (size_t g=0; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    size_t cnt = 0;
                    for (size_t l=0; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        cnt += (fabs(*X1-*X2)>thresh);
                    }
                    *Y = cnt;
                }
            }
        }
    }

    return 0;
}


int dist0_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in dist0_c: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist0_c: vectors in X1 and X2 must have the same length\n"); return 1; }
    float thresh = 2.0f * FLT_EPSILON;

    if (N==0) {}
    else if (L==N)
    {
        size_t cnt = 0;
        for (size_t l=0; l<L; ++l)
        {
            cnt += (fabsf(*X1++-*X2++)>thresh || fabsf(*X1++-*X2++)>thresh);
        }
        *Y = cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? 2*L : 0, J2 = (L==N2) ? 2*L : 0;
            for (size_t v=0; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                size_t cnt = 0;
                for (size_t l=0; l<L; ++l)
                {
                    cnt += (fabsf(*X1++-*X2++)>thresh || fabsf(*X1++-*X2++)>thresh);
                }
                *Y = cnt;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 2, J2 = (L==N2) ? 0 : 2;
            const size_t K1 = (L==N1) ? 2 : 2*K, K2 = (L==N2) ? 2 : 2*K;
            const size_t I1 = (L==N1) ? 0 : 2*B*(L-1), I2 = (L==N2) ? 0 : 2*B*(L-1);
            for (size_t g=0; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    size_t cnt = 0;
                    for (size_t l=0; l<L; ++l, X1+=K1-2, X2+=K2-2)
                    {
                        cnt += (fabsf(*X1++-*X2++)>thresh || fabsf(*X1++-*X2++)>thresh);
                    }
                    *Y = cnt;
                }
            }
        }
    }

    return 0;
}


int dist0_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in dist0_z: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist0_z: vectors in X1 and X2 must have the same length\n"); return 1; }
    double thresh = 2.0 * DBL_EPSILON;

    if (N==0) {}
    else if (L==N)
    {
        size_t cnt = 0;
        for (size_t l=0; l<L; ++l)
        {
            cnt += (fabs(*X1++-*X2++)>thresh || fabs(*X1++-*X2++)>thresh);
        }
        *Y = cnt;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            const size_t J1 = (L==N1) ? 2*L : 0, J2 = (L==N2) ? 2*L : 0;
            for (size_t v=0; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                size_t cnt = 0;
                for (size_t l=0; l<L; ++l)
                {
                    cnt += (fabs(*X1++-*X2++)>thresh || fabs(*X1++-*X2++)>thresh);
                }
                *Y = cnt;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0 : 2, J2 = (L==N2) ? 0 : 2;
            const size_t K1 = (L==N1) ? 2 : 2*K, K2 = (L==N2) ? 2 : 2*K;
            const size_t I1 = (L==N1) ? 0 : 2*B*(L-1), I2 = (L==N2) ? 0 : 2*B*(L-1);
            for (size_t g=0; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    size_t cnt = 0;
                    for (size_t l=0; l<L; ++l, X1+=K1-2, X2+=K2-2)
                    {
                        cnt += (fabs(*X1++-*X2++)>thresh || fabs(*X1++-*X2++)>thresh);
                    }
                    *Y = cnt;
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
