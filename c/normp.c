//Vec2scalar (reduction) operation.
//Gets the p-norm for each vector in X along dim.
//This is the Lp norm of each vector in X.
//For each vector, y = sum(|x|^p)^1/p.
//For complex case, output is real.

#include <stdio.h>
#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int normp_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p)
{
    if (dim>3u) { fprintf(stderr,"error in normp_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float ip = 1.0f / p;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        float sm = 0.0f;
        for (size_t l=L; l>0u; --l, ++X) { sm += powf(fabsf(*X),p); }
        *Y = powf(sm,ip);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float sm;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                sm = 0.0f;
                for (size_t l=L; l>0u; --l, ++X) { sm += powf(fabsf(*X),p); }
                *Y = powf(sm,ip);
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = powf(fabsf(*X),p); }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += powf(fabsf(*X),p); }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y = powf(*Y,ip); }
        }
        else
        {
            float sm;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=K) { sm += powf(fabsf(*X),p); }
                    *Y = powf(sm,ip);
                }
            }
        }
    }

    return 0;
}


int normp_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p)
{
    if (dim>3u) { fprintf(stderr,"error in normp_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double ip = 1.0 / p;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        double sm = 0.0;
        for (size_t l=L; l>0u; --l, ++X) { sm += pow(fabs(*X),p); }
        *Y = pow(sm,ip);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double sm;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                sm = 0.0;
                for (size_t l=L; l>0u; --l, ++X) { sm += pow(fabs(*X),p); }
                *Y = pow(sm,ip);
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = pow(fabs(*X),p); }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += pow(fabs(*X),p); }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y = pow(*Y,ip); }
        }
        else
        {
            double sm;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=L; l>0u; --l, X+=K) { sm += pow(fabs(*X),p); }
                    *Y = pow(sm,ip);
                }
            }
        }
    }

    return 0;
}


int normp_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const float p)
{
    if (dim>3u) { fprintf(stderr,"error in normp_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float ip = 1.0f / p;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, X+=2, ++Y) { *Y = sqrtf(*X**X + *(X+1)**(X+1)); }
    }
    else if (L==N)
    {
        float sm = 0.0f;
        for (size_t l=L; l>0u; --l, X+=2) { sm += powf(sqrtf(*X**X + *(X+1)**(X+1)),p); }
        *Y = powf(sm,ip);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float sm;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                sm = 0.0f;
                for (size_t l=L; l>0u; --l, X+=2) { sm += powf(sqrtf(*X**X + *(X+1)**(X+1)),p); }
                *Y = powf(sm,ip);
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y = powf(sqrtf(*X**X + *(X+1)**(X+1)),p); }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y += powf(sqrtf(*X**X + *(X+1)**(X+1)),p); }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y = powf(*Y,ip); }
        }
        else
        {
            float sm;
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=2u*K) { sm += powf(sqrtf(*X**X + *(X+1)**(X+1)),p); }
                    *Y = powf(sm,ip);
                }
            }
        }
    }

    return 0;
}


int normp_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const double p)
{
    if (dim>3u) { fprintf(stderr,"error in normp_z: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double ip = 1.0 / p;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, X+=2, ++Y) { *Y = sqrt(*X**X + *(X+1)**(X+1)); }
    }
    else if (L==N)
    {
        double sm = 0.0;
        for (size_t l=L; l>0u; --l, X+=2) { sm += pow(sqrt(*X**X + *(X+1)**(X+1)),p); }
        *Y = pow(sm,ip);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double sm;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                sm = 0.0;
                for (size_t l=L; l>0u; --l, X+=2) { sm += pow(sqrt(*X**X + *(X+1)**(X+1)),p); }
                *Y = pow(sm,ip);
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y = pow(sqrt(*X**X + *(X+1)**(X+1)),p); }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, X+=2, ++Y) { *Y += pow(sqrt(*X**X + *(X+1)**(X+1)),p); }
            }
            for (size_t v=V; v>0u; --v, ++Y) { *Y = pow(*Y,ip); }
        }
        else
        {
            double sm;
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=L; l>0u; --l, X+=2u*K) { sm += pow(sqrt(*X**X + *(X+1)**(X+1)),p); }
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
