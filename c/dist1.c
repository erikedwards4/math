//Vecs2scalar operation for 2 inputs X1 and X2.
//Gets L1 distance between pairs of vecs in X1, X2.

//This is the 'Manhattan' or 'taxicab' or 'rectlinear' vector distance,
//which is just the sum of absolute differences: sum(|X1-X2|).

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int dist1_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int dist1_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int dist1_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int dist1_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);


int dist1_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dist1_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist1_s: vectors in X1 and X2 must have the same length\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X1, ++X2, ++Y) { *Y = fabsf(*X1-*X2); }
    }
    else if (L==N)
    {
        float sm = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            sm += fabsf(*X1-*X2);
        }
        *Y = sm;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? L : 0u, J2 = (L==N2) ? L : 0u;
            float sm;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                sm = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    sm += fabsf(*X1-*X2);
                }
                *Y = sm;
            }
        }
        else if (G==1u)
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, ++Y) { *Y = fabsf(*X1-*X2); }
            Y -= V; X1 += 1u-J1; X2 += 1u-J2;
            for (size_t l=1u; l<L; ++l, Y-=V, X1+=1u-J1, X2+=1u-J2)
            {
                for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, ++Y) { *Y += fabsf(*X1-*X2); }
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : B*(L-1u), I2 = (L==N2) ? 0u : B*(L-1u);
            float sm;
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        sm += fabsf(*X1-*X2);
                    }
                    *Y = sm;
                }
            }
        }
    }

    return 0;
}


int dist1_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dist1_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist1_s: vectors in X1 and X2 must have the same length\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X1, ++X2, ++Y) { *Y = fabs(*X1-*X2); }
    }
    else if (L==N)
    {
        double sm = 0.0;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            sm += fabs(*X1-*X2);
        }
        *Y = sm;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? L : 0u, J2 = (L==N2) ? L : 0u;
            double sm;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                sm = 0.0;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    sm += fabs(*X1-*X2);
                }
                *Y = sm;
            }
        }
        else if (G==1u)
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, ++Y) { *Y = fabs(*X1-*X2); }
            Y -= V; X1 += 1u-J1; X2 += 1u-J2;
            for (size_t l=1u; l<L; ++l, Y-=V, X1+=1u-J1, X2+=1u-J2)
            {
                for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, ++Y) { *Y += fabs(*X1-*X2); }
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : B*(L-1u), I2 = (L==N2) ? 0u : B*(L-1u);
            double sm;
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2)
                    {
                        sm += fabs(*X1-*X2);
                    }
                    *Y = sm;
                }
            }
        }
    }

    return 0;
}


int dist1_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dist1_c: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist1_c: vectors in X1 and X2 must have the same length\n"); return 1; }
    float dr, di;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X1, ++X2, ++Y)
        {
            dr = *X1 - *X2;
            di = *++X1 - *++X2;
            *Y = sqrtf(dr*dr+di*di);
        }
    }
    else if (L==N)
    {
        float sm = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            dr = *X1 - *X2;
            di = *++X1 - *++X2;
            sm += sqrtf(dr*dr+di*di);
        }
        *Y = sm;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? 2u*L : 0u, J2 = (L==N2) ? 2u*L : 0u;
            float sm;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                sm = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    dr = *X1 - *X2;
                    di = *++X1 - *++X2;
                    sm += sqrtf(dr*dr+di*di);
                }
                *Y = sm;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 2u, J2 = (L==N2) ? 0u : 2u;
            const size_t K1 = (L==N1) ? 2u : 2u*K, K2 = (L==N2) ? 2u : 2u*K;
            const size_t I1 = (L==N1) ? 0u : 2u*B*(L-1u), I2 = (L==N2) ? 0u : 2u*B*(L-1u);
            float sm;
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0u; l<L; ++l, X1+=K1-1u, X2+=K2-1u)
                    {
                        dr = *X1 - *X2;
                        di = *++X1 - *++X2;
                        sm += sqrtf(dr*dr+di*di);
                    }
                    *Y = sm;
                }
            }
        }
    }

    return 0;
}


int dist1_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dist1_z: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2*C2*S2*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dist1_z: vectors in X1 and X2 must have the same length\n"); return 1; }
    double dr, di;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X1, ++X2, ++Y)
        {
            dr = *X1 - *X2;
            di = *++X1 - *++X2;
            *Y = sqrt(dr*dr+di*di);
        }
    }
    else if (L==N)
    {
        double sm = 0.0;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2)
        {
            dr = *X1 - *X2;
            di = *++X1 - *++X2;
            sm += sqrt(dr*dr+di*di);
        }
        *Y = sm;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? 2u*L : 0u, J2 = (L==N2) ? 2u*L : 0u;
            double sm;
            for (size_t v=0u; v<V; ++v, X1-=J1, X2-=J2, ++Y)
            {
                sm = 0.0;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2)
                {
                    dr = *X1 - *X2;
                    di = *++X1 - *++X2;
                    sm += sqrt(dr*dr+di*di);
                }
                *Y = sm;
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 2u, J2 = (L==N2) ? 0u : 2u;
            const size_t K1 = (L==N1) ? 2u : 2u*K, K2 = (L==N2) ? 2u : 2u*K;
            const size_t I1 = (L==N1) ? 0u : 2u*B*(L-1u), I2 = (L==N2) ? 0u : 2u*B*(L-1u);
            double sm;
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=L*K1-J1, X2-=L*K2-J2, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0u; l<L; ++l, X1+=K1-1u, X2+=K2-1u)
                    {
                        dr = *X1 - *X2;
                        di = *++X1 - *++X2;
                        sm += sqrt(dr*dr+di*di);
                    }
                    *Y = sm;
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
