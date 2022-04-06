//Vec2scalar (reduction) operation.
//Gets product of elements for each vector in X along dim.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int prod_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in prod_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        float pr = *X++;
        for (size_t l=L; l>1u; --l, ++X) { pr *= *X; }
        *Y = pr;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float pr;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                pr = *X++;
                for (size_t l=L; l>1u; --l, ++X) { pr *= *X; }
                *Y = pr;
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y *= *X; }
            }
        }
        else
        {
            float pr;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    pr = *X; X += K;
                    for (size_t l=L; l>1u; --l, X+=K) { pr *= *X; }
                    *Y = pr;
                }
            }
        }
    }

    return 0;
}


int prod_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in prod_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        double pr = *X++;
        for (size_t l=L; l>1u; --l, ++X) { pr *= *X; }
        *Y = pr;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double pr;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                pr = *X++;
                for (size_t l=L; l>1u; --l, ++X) { pr *= *X; }
                *Y = pr;
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=L; l>1u; --l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y *= *X; }
            }
        }
        else
        {
            double pr;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    pr = *X; X += K;
                    for (size_t l=L; l>1u; --l, X+=K) { pr *= *X; }
                    *Y = pr;
                }
            }
        }
    }

    return 0;
}


int prod_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in prod_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float yr, yi;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y = *X++; *(Y+1) = *X++;
        for (size_t l=L; l>1u; --l, X+=2)
        {
            yr = *X**Y - *(X+1)**(Y+1);
            yi = *X**(Y+1) + *(X+1)**Y;
            *Y = yr; *(Y+1) = yi;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, Y+=2)
            {
                *Y = *X++; *(Y+1) = *X++;
                for (size_t l=L; l>1u; --l, X+=2)
                {
                    yr = *X**Y - *(X+1)**(Y+1);
                    yi = *X**(Y+1) + *(X+1)**Y;
                    *Y = yr; *(Y+1) = yi;
                }
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            Y -= 2u*V;
            for (size_t l=L; l>1u; --l, Y-=2u*V)
            {
                for (size_t v=V; v>0u; --v, X+=2)
                {
                    yr = *X**Y - *(X+1)**(Y+1);
                    yi = *X**(Y+1) + *(X+1)**Y;
                    *Y++ = yr; *Y++ = yi;
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y+=2)
                {
                    *Y = *X; *(Y+1) = *(X+1); X += 2u*K;
                    for (size_t l=L; l>1u; --l, X+=2u*K)
                    {
                        yr = *X**Y - *(X+1)**(Y+1);
                        yi = *X**(Y+1) + *(X+1)**Y;
                        *Y = yr; *(Y+1) = yi;
                    }
                }
            }
        }
    }
    
    return 0;
}


int prod_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in prod_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double yr, yi;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y = *X++; *(Y+1) = *X++;
        for (size_t l=L; l>1u; --l, X+=2)
        {
            yr = *X**Y - *(X+1)**(Y+1);
            yi = *X**(Y+1) + *(X+1)**Y;
            *Y = yr; *(Y+1) = yi;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, Y+=2)
            {
                *Y = *X++; *(Y+1) = *X++;
                for (size_t l=L; l>1u; --l, X+=2)
                {
                    yr = *X**Y - *(X+1)**(Y+1);
                    yi = *X**(Y+1) + *(X+1)**Y;
                    *Y = yr; *(Y+1) = yi;
                }
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            Y -= 2u*V;
            for (size_t l=L; l>1u; --l, Y-=2u*V)
            {
                for (size_t v=V; v>0u; --v, X+=2)
                {
                    yr = *X**Y - *(X+1)**(Y+1);
                    yi = *X**(Y+1) + *(X+1)**Y;
                    *Y++ = yr; *Y++ = yi;
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y+=2)
                {
                    *Y = *X; *(Y+1) = *(X+1); X += 2u*K;
                    for (size_t l=L; l>1u; --l, X+=2u*K)
                    {
                        yr = *X**Y - *(X+1)**(Y+1);
                        yi = *X**(Y+1) + *(X+1)**Y;
                        *Y = yr; *(Y+1) = yi;
                    }
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
