//Vecs2scalar operation for 2 inputs X1 and X2.
//Dot product for each pair of vectors.

//For complex inputs, this uses the unconjugated dot product.

//Direct for loop is much faster for small N, faster for medium N,
//and ~same speed for huge N (~1e6), compared to cblas_?dot.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int dot_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int dot_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int dot_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);
int dot_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim);


int dot_s (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dot_s: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2u*C2u*S2u*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dot_s: vectors in X1 and X2 must have the same length\n"); return 1; }

    if (N==0u) {}
    else if (L==N)
    {
        float sm2 = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2) { sm2 += *X1 * *X2; }
        *Y = sm2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float sm2;
            const int J1 = (L==N1) ? -(int)L : 0, J2 = (L==N2) ? -(int)L : 0;
            for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, ++Y)
            {
                sm2 = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2) { sm2 += *X1 * *X2; }
                *Y = sm2;
            }
        }
        else if (G==1u)
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, ++Y) { *Y = *X1 * *X2; }
            X1 += 1u-J1; X2 += 1u-J2;  Y -= V;
            for (size_t l=1u; l<L; ++l, X1+=1u-J1, X2+=1u-J2, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, ++Y) { *Y += *X1 * *X2; }
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : B*(L-1u), I2 = (L==N2) ? 0u : B*(L-1u);
            float sm2;
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=K1*L-J1, X2-=K2u*L-J2, ++Y)
                {
                    sm2 = 0.0f;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2) { sm2 += *X1 * *X2; }
                    *Y = sm2;
                }
            }
        }
    }

    return 0;
}


int dot_d (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dot_d: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2u*C2u*S2u*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dot_d: vectors in X1 and X2 must have the same length\n"); return 1; }

    if (N==0u) {}
    else if (L==N)
    {
        double sm2 = 0.0;
        for (size_t l=0u; l<L; ++l, ++X1, ++X2) { sm2 += *X1 * *X2; }
        *Y = sm2;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double sm2;
            const int J1 = (L==N1) ? -(int)L : 0, J2 = (L==N2) ? -(int)L : 0;
            for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, ++Y)
            {
                sm2 = 0.0;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2) { sm2 += *X1 * *X2; }
                *Y = sm2;
            }
        }
        else if (G==1u)
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, ++Y) { *Y = *X1 * *X2; }
            X1 += 1u-J1; X2 += 1u-J2;  Y -= V;
            for (size_t l=1u; l<L; ++l, X1+=1u-J1, X2+=1u-J2, Y-=V)
            {
                for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, ++Y) { *Y += *X1 * *X2; }
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 1u, J2 = (L==N2) ? 0u : 1u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : B*(L-1u), I2 = (L==N2) ? 0u : B*(L-1u);
            double sm2;
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1-=K1*L-J1, X2-=K2u*L-J2, ++Y)
                {
                    sm2 = 0.0;
                    for (size_t l=0u; l<L; ++l, X1+=K1, X2+=K2) { sm2 += *X1 * *X2; }
                    *Y = sm2;
                }
            }
        }
    }

    return 0;
}


int dot_c (float *Y, const float *X1, const float *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dot_c: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2u*C2u*S2u*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dot_c: vectors in X1 and X2 must have the same length\n"); return 1; }

    if (N==0u) {}
    else if (L==N)
    {
        if (L<5000u)
        {
            float sm2r = 0.0f, sm2i = 0.0f;
            for (size_t l=0u; l<L; ++l, X1+=2, X2+=2) { sm2r += *X1**X2 - *(X1+1)**(X2+1); sm2i += *X1**(X2+1) + *(X1+1)**X2; }
            *Y = sm2r; *++Y = sm2i;
        }
        else { cblas_cdotu_sub((int)L,X1,1,X2,1,(_Complex float *)Y); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? 0u : 2u*L, J2 = (L==N2) ? 0u : 2u*L;
            for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, Y+=2)
            {
                cblas_cdotu_sub((int)L,X1,1,X2,1,(_Complex float *)Y);
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 2u, J2 = (L==N2) ? 0u : 2u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : 2u*B*(L-1u), I2 = (L==N2) ? 0u : 2u*B*(L-1u);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1+=J1, X2+=J2, Y+=2)
                {
                    cblas_cdotu_sub((int)L,X1,(int)K1,X2,(int)K2,(_Complex float *)Y);
                }
            }
        }
    }

    return 0;
}


int dot_z (double *Y, const double *X1, const double *X2, const size_t R1, const size_t C1, const size_t S1, const size_t H1, const size_t R2, const size_t C2, const size_t S2, const size_t H2, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in dot_z: dim must be in [0 3]\n"); return 1; }

    const size_t R = (R1>R2) ? R1 : R2;
    const size_t C = (C1>C2) ? C1 : C2;
    const size_t S = (S1>S2) ? S1 : S2;
    const size_t H = (H1>H2) ? H1 : H2;
    const size_t N = R*C*S*H, N1 = R1*C1*S1*H1, N2 = R2u*C2u*S2u*H2;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L1 = (dim==0u) ? R1 : (dim==1u) ? C1 : (dim==2u) ? S1 : H1;
    const size_t L2 = (dim==0u) ? R2 : (dim==1u) ? C2 : (dim==2u) ? S2 : H2;
    if (L1!=L2) { fprintf(stderr,"error in dot_z: vectors in X1 and X2 must have the same length\n"); return 1; }

    if (N==0u) {}
    else if (L==N)
    {
        if (L<5000u)
        {
            double sm2r = 0.0, sm2i = 0.0;
            for (size_t l=0u; l<L; ++l, X1+=2, X2+=2) { sm2r += *X1**X2 - *(X1+1)**(X2+1); sm2i += *X1**(X2+1) + *(X1+1)**X2; }
            *Y = sm2r; *++Y = sm2i;
        }
        else { cblas_zdotu_sub((int)L,X1,1,X2,1,(_Complex double *)Y); }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            const size_t J1 = (L==N1) ? 0u : 2u*L, J2 = (L==N2) ? 0u : 2u*L;
            for (size_t v=0u; v<V; ++v, X1+=J1, X2+=J2, Y+=2)
            {
                cblas_zdotu_sub((int)L,X1,1,X2,1,(_Complex double *)Y);
            }
        }
        else
        {
            const size_t J1 = (L==N1) ? 0u : 2u, J2 = (L==N2) ? 0u : 2u;
            const size_t K1 = (L==N1) ? 1u : K, K2 = (L==N2) ? 1u : K;
            const size_t I1 = (L==N1) ? 0u : 2u*B*(L-1u), I2 = (L==N2) ? 0u : 2u*B*(L-1u);
            for (size_t g=0u; g<G; ++g, X1+=I1, X2+=I2)
            {
                for (size_t b=0u; b<B; ++b, X1+=J1, X2+=J2, Y+=2)
                {
                    cblas_zdotu_sub((int)L,X1,(int)K1,X2,(int)K2,(_Complex double *)Y);
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
