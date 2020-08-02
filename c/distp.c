//Vecs2scalar operation for 2 inputs X1 and X2.
//Gets Lp distance between pairs of vecs in X1, X2.

//This sums the pth power of the difference vecs,
//and then takes the pth root: sum(|X1-X2|^p)^1/p.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int distp_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const float p);
int distp_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const double p);
int distp_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const float p);
int distp_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const double p);


int distp_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const float p)
{
    if (dim>3) { fprintf(stderr,"error in distp_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in distp_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float ip = 1.0f / p;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X1, ++X2, ++Y) { *Y = fabsf(*X1-*X2); }
    }
    else if (L==N)
    {
        float sm = 0.0f;
        for (size_t l=0; l<L; ++l, ++X1, ++X2)
        {
            sm += powf(fabsf(*X1-*X2),p);
        }
        *Y = powf(sm,ip);
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
                float sm = 0.0f;
                for (size_t l=0; l<L; ++l, ++X1, ++X2)
                {
                    sm += powf(fabsf(*X1-*X2),p);
                }
                *Y = powf(sm,ip);
            }
        }
        else if (G==1)
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            for (size_t v=0; v<V; ++v, X1+=J1, X2+=J2, ++Y) { *Y = powf(fabsf(*X1-*X2),p); }
            Y -= V; X1 += 1-J1; X2 += 1-J2;
            for (size_t l=1; l<L; ++l, Y-=V, X1+=1-J1, X2+=1-J2)
            {
                for (size_t v=0; v<V; ++v, X1+=J1, X2+=J2, ++Y) { *Y += powf(fabsf(*X1-*X2),p); }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = powf(*Y,ip); }
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
                    float sm = 0.0f;
                    for (size_t l=0; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        sm += powf(fabsf(*X1-*X2),p);
                    }
                    *Y = powf(sm,ip);
                }
            }
        }
    }

    return 0;
}


int distp_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const double p)
{
    if (dim>3) { fprintf(stderr,"error in distp_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in distp_s: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double ip = 1.0 / p;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X1, ++X2, ++Y) { *Y = fabs(*X1-*X2); }
    }
    else if (L==N)
    {
        double sm = 0.0;
        for (size_t l=0; l<L; ++l, ++X1, ++X2)
        {
            sm += pow(fabs(*X1-*X2),p);
        }
        *Y = pow(sm,ip);
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
                double sm = 0.0;
                for (size_t l=0; l<L; ++l, ++X1, ++X2)
                {
                    sm += pow(fabs(*X1-*X2),p);
                }
                *Y = pow(sm,ip);
            }
        }
        else if (G==1)
        {
            const size_t J1 = (L==N1) ? 0 : 1, J2 = (L==N2) ? 0 : 1;
            for (size_t v=0; v<V; ++v, X1+=J1, X2+=J2, ++Y) { *Y = pow(fabs(*X1-*X2),p); }
            Y -= V; X1 += 1-J1; X2 += 1-J2;
            for (size_t l=1; l<L; ++l, Y-=V, X1+=1-J1, X2+=1-J2)
            {
                for (size_t v=0; v<V; ++v, X1+=J1, X2+=J2, ++Y) { *Y += pow(fabs(*X1-*X2),p); }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = pow(*Y,ip); }
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
                    double sm = 0.0;
                    for (size_t l=0; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        sm += pow(fabs(*X1-*X2),p);
                    }
                    *Y = pow(sm,ip);
                }
            }
        }
    }

    return 0;
}


int distp_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const float p)
{
    if (dim>3) { fprintf(stderr,"error in distp_c: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in distp_c: vectors in X1 and X2 must have the same length\n"); return 1; }
    const float p2 = p/2.0f, ip = 1.0f/p;
    float dr, di;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++Y)
        {
            dr = *X1++ - *X2++;
            di = *X1++ - *X2++;
            *Y = sqrtf(dr*dr+di*di);
        }
    }
    else if (L==N)
    {
        float sm = 0.0f;
        for (size_t l=0; l<L; ++l)
        {
            dr = *X1++ - *X2++;
            di = *X1++ - *X2++;
            sm += powf(dr*dr+di*di,p2);
        }
        *Y = powf(sm,ip);
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
                float sm = 0.0f;
                for (size_t l=0; l<L; ++l)
                {
                    dr = *X1++ - *X2++;
                    di = *X1++ - *X2++;
                    sm += powf(dr*dr+di*di,p2);
                }
                *Y = powf(sm,ip);
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
                    float sm = 0.0;
                    for (size_t l=0; l<L; ++l, X1+=K1-1, X2+=K2-1)
                    {
                        dr = *X1++ - *X2++;
                        di = *X1 - *X2;
                        sm += powf(dr*dr+di*di,p2);
                    }
                    *Y = powf(sm,ip);
                }
            }
        }
    }

    return 0;
}


int distp_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim, const double p)
{
    if (dim>3) { fprintf(stderr,"error in distp_z: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t L1 = (dim==0) ? R1 : (dim==1) ? C1 : (dim==2) ? S1 : H1;
    const size_t L2 = (dim==0) ? R2 : (dim==1) ? C2 : (dim==2) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in distp_z: vectors in X1 and X2 must have the same length\n"); return 1; }
    const double p2 = p/2.0, ip = 1.0/p;
    double dr, di;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++Y)
        {
            dr = *X1++ - *X2++;
            di = *X1++ - *X2++;
            *Y = sqrt(dr*dr+di*di);
        }
    }
    else if (L==N)
    {
        double sm = 0.0;
        for (size_t l=0; l<L; ++l)
        {
            dr = *X1++ - *X2++;
            di = *X1++ - *X2++;
            sm += pow(dr*dr+di*di,p2);
        }
        *Y = pow(sm,ip);
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
                double sm = 0.0;
                for (size_t l=0; l<L; ++l)
                {
                    dr = *X1++ - *X2++;
                    di = *X1++ - *X2++;
                    sm += pow(dr*dr+di*di,p2);
                }
                *Y = pow(sm,ip);
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
                    double sm = 0.0;
                    for (size_t l=0; l<L; ++l, X1+=K1-1, X2+=K2-1)
                    {
                        dr = *X1++ - *X2++;
                        di = *X1 - *X2;
                        sm += pow(dr*dr+di*di,p2);
                    }
                    *Y = pow(sm,ip);
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
