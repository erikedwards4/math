//Vec2scalar (reduction) operation.
//Gets the coefficient of variation (std/mean) for each vector in X along dim.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int coeff_var_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);
int coeff_var_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased);


int coeff_var_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased)
{
    if (dim>3u) { fprintf(stderr,"error in coeff_var_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float den = 1.0f/(float)L, den2 = (biased) ? den : 1.0f/(float)(L-1u);

    if (N==0u) {}
    else if (L<2u) { fprintf(stderr,"error in coeff_var_s: L must be > 1\n"); return 1; }
    else if (L==N)
    {
        float x, mn = 0.0f, sm2 = 0.0f;
        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=0u; l<L; ++l) { x = *--X - mn; sm2 += x*x; }
        *Y = sqrtf(sm2*den2) / mn;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float x, mn, sm2;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                mn = sm2 = 0.0f;
                for (size_t l=0u; l<L; ++l, ++X) { *Y += *X; }
                mn *= den;
                X -= L;
                for (size_t l=0u; l<L; ++l, ++X) { x = *X - mn; sm2 += x*x; }
                *Y = sqrtf(sm2*den2) / mn;
            }
        }
        else if (G==1u)
        {
            float x, *mn;
            if (!(mn=(float *)calloc(V,sizeof(float)))) { fprintf(stderr,"error in coeff_var_s: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0u; l<L; ++l, mn-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++mn) { *mn += *X; }
            }
            X -= N;
            for (size_t v=V; v>0u; --v, ++mn, ++Y) { *mn *= den; *Y = 0.0f; }
            mn -= V; Y -= V;
            for (size_t l=0u; l<L; ++l, mn-=V, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++mn, ++Y) { x = *X - *mn; *Y += x*x; }
            }
            for (size_t v=V; v>0u; --v, ++mn, ++Y) { *Y = sqrtf(*Y*den2) / *mn; }
            mn -= V; free(mn);
        }
        else
        {
            float x, mn, sm2;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X, ++Y)
                {
                    mn = sm2 = 0.0f;
                    for (size_t l=0u; l<L-1u; ++l, X+=K) { mn += *X; }
                    mn += *X; mn *= den;
                    for (size_t l=0u; l<L-1u; ++l, X-=K) { x = *X - mn; sm2 += x*x; }
                    x = *X - mn; sm2 += x*x;
                    *Y = sqrtf(sm2*den2) / mn;
                }
            }
        }
    }

    return 0;
}


int coeff_var_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int biased)
{
    if (dim>3u) { fprintf(stderr,"error in coeff_var_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double den = 1.0/(double)L, den2 = (biased) ? den : 1.0/(double)(L-1u);

    if (N==0u) {}
    else if (L<2u) { fprintf(stderr,"error in coeff_var_d: L must be > 1\n"); return 1; }
    else if (L==N)
    {
        double x, mn = 0.0, sm2 = 0.0;
        for (size_t l=0u; l<L; ++l, ++X) { mn += *X; }
        mn *= den;
        for (size_t l=0u; l<L; ++l) { x = *--X - mn; sm2 += x*x; }
        *Y = sqrt(sm2*den2) / mn;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double x, mn, sm2;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                mn = sm2 = 0.0;
                for (size_t l=0u; l<L; ++l, ++X) { *Y += *X; }
                mn *= den;
                X -= L;
                for (size_t l=0u; l<L; ++l, ++X) { x = *X - mn; sm2 += x*x; }
                *Y = sqrt(sm2*den2) / mn;
            }
        }
        else if (G==1u)
        {
            double x, *mn;
            if (!(mn=(double *)calloc(V,sizeof(double)))) { fprintf(stderr,"error in coeff_var_d: problem with calloc. "); perror("calloc"); return 1; }
            for (size_t l=0u; l<L; ++l, mn-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++mn) { *mn += *X; }
            }
            X -= N;
            for (size_t v=V; v>0u; --v, ++mn, ++Y) { *mn *= den; *Y = 0.0; }
            mn -= V; Y -= V;
            for (size_t l=0u; l<L; ++l, mn-=V, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++mn, ++Y) { x = *X - *mn; *Y += x*x; }
            }
            for (size_t v=V; v>0u; --v, ++mn, ++Y) { *Y = sqrt(*Y*den2) / *mn; }
            mn -= V; free(mn);
        }
        else
        {
            double x, mn, sm2;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, ++X, ++Y)
                {
                    mn = sm2 = 0.0;
                    for (size_t l=0u; l<L-1u; ++l, X+=K) { mn += *X; }
                    mn += *X; mn *= den;
                    for (size_t l=0u; l<L-1u; ++l, X-=K) { x = *X - mn; sm2 += x*x; }
                    x = *X - mn; sm2 += x*x;
                    *Y = sqrt(sm2*den2) / mn;
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
