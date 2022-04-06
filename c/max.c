//Vec2scalar (reduction) operation.
//Gets maximum of values for each vector in X along dim.
//For complex case, finds max absolute value and outputs the complex number.
//For complex case, this is the proper absolute value; see amax for the other one.

//Computational note: it is much faster to use intermediate mx variable (on stack) than to compare to *Y.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int max_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in max_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float mx;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        mx = *X++;
        for (size_t l=1u; l<L; ++l, ++X) { if (*X>mx) { mx = *X; } }
        *Y = mx;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++Y)
            {
                mx = *X++;
                for (size_t l=1u; l<L; ++l, ++X) { if (*X>mx) { mx = *X; } }
                *Y = mx;
            }
        }
        // else if (G==1u)
        // {
        //     for (size_t l=L; l>0u; --l, Y-=V)
        //     {
        //         for (size_t v=V; v>0u; --v, ++X, ++Y)
        //         {
        //             if (l==0u || *X>*Y) { *Y = *X; }
        //         }
        //     }
        // }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    mx = *X; X += K;
                    for (size_t l=1u; l<L; ++l, X+=K) { if (*X>mx) { mx = *X; } }
                    *Y = mx;
                }
            }
        }
    }

    return 0;
}


int max_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in max_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double mx;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        mx = *X++;
        for (size_t l=1u; l<L; ++l, ++X) { if (*X>mx) { mx = *X; } }
        *Y = mx;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++Y)
            {
                mx = *X++;
                for (size_t l=1u; l<L; ++l, ++X) { if (*X>mx) { mx = *X; } }
                *Y = mx;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    mx = *X; X += K;
                    for (size_t l=1u; l<L; ++l, X+=K) { if (*X>mx) { mx = *X; } }
                    *Y = mx;
                }
            }
        }
    }

    return 0;
}


int max_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in max_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float xx, mx;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        mx = *X**X + *(X+1)**(X+1);
        *Y = *X; *(Y+1) = *(X+1); X += 2;
        for (size_t l=1u; l<L; ++l, X+=2)
        {
            xx = *X**X + *(X+1)**(X+1);
            if (xx>mx) { mx = xx; *Y = *X; *(Y+1) = *(X+1); }
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
                mx = *X**X + *(X+1)**(X+1);
                *Y = *X; *(Y+1) = *(X+1); X += 2;
                for (size_t l=1u; l<L; ++l, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx>mx) { mx = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y+=2)
                {
                    mx = *X**X + *(X+1)**(X+1);
                    *Y = *X; *(Y+1) = *(X+1); X += 2u*K;
                    for (size_t l=1u; l<L; ++l, X+=2u*K)
                    {
                        xx = *X**X + *(X+1)**(X+1);
                        if (xx>mx) { mx = xx; *Y = *X; *(Y+1) = *(X+1); }
                    }
                }
            }
        }
    }

    return 0;
}


int max_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in max_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double xx, mx;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        mx = *X**X + *(X+1)**(X+1);
        *Y = *X; *(Y+1) = *(X+1); X += 2;
        for (size_t l=1u; l<L; ++l, X+=2)
        {
            xx = *X**X + *(X+1)**(X+1);
            if (xx>mx) { mx = xx; *Y = *X; *(Y+1) = *(X+1); }
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
                mx = *X**X + *(X+1)**(X+1);
                *Y = *X; *(Y+1) = *(X+1); X += 2;
                for (size_t l=1u; l<L; ++l, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx>mx) { mx = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y+=2)
                {
                    mx = *X**X + *(X+1)**(X+1);
                    *Y = *X; *(Y+1) = *(X+1); X += 2u*K;
                    for (size_t l=1u; l<L; ++l, X+=2u*K)
                    {
                        xx = *X**X + *(X+1)**(X+1);
                        if (xx>mx) { mx = xx; *Y = *X; *(Y+1) = *(X+1); }
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
