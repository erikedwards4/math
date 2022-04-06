//Flips X along dim (reverses order of elements).
//For d=0, this is flipud. For d=1, this is fliplr.
//This has in-place and not-in-place versions.

//This is a vec2vec operation (see readme.txt).

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int flip_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in flip_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        Y += L - 1u;
        for (size_t l=L; l>0u; --l, ++X, --Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            Y += L - 1u;
            for (size_t v=V; v>0u; --v, Y+=2u*L)
            {
                for (size_t l=L; l>0u; --l, ++X, --Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            Y += V*(L-1u);
            for (size_t l=L; l>0u; --l, Y-=2u*V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y+=K+1u)
                {
                    Y += K*(L-1u);
                    for (size_t l=L; l>0u; --l, Y-=K, X+=K) { *Y = *X; }
                }
            }
        }
    }

    return 0;
}


int flip_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in flip_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        Y += L - 1u;
        for (size_t l=L; l>0u; --l, ++X, --Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            Y += L - 1u;
            for (size_t v=V; v>0u; --v, Y+=2u*L)
            {
                for (size_t l=L; l>0u; --l, ++X, --Y) { *Y = *X; }
            }
        }
        else if (G==1u)
        {
            Y += V*(L-1u);
            for (size_t l=L; l>0u; --l, Y-=2u*V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, Y+=K+1u)
                {
                    Y += K*(L-1u);
                    for (size_t l=L; l>0u; --l, Y-=K, X+=K) { *Y = *X; }
                }
            }
        }
    }

    return 0;
}


int flip_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in flip_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        Y += 2u*(L-1u);
        for (size_t l=L; l>0u; --l, ++X, Y-=2) { *Y = *X; *(Y+1) = *++X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            Y += 2u*(L-1u);
            for (size_t v=V; v>0u; --v, Y+=4u*L)
            {
                for (size_t l=L; l>0u; --l, ++X, Y-=2) { *Y = *X; *(Y+1) = *++X; }
            }
        }
        else if (G==1u)
        {
            Y += 2u*V*(L-1u);
            for (size_t l=L; l>0u; --l, Y-=4u*V)
            {
                for (size_t v=0u; v<2u*V; ++v, ++X, ++Y) { *Y = *X; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y+=2u*K+2u)
                {
                    Y += 2u*K*(L-1u);
                    for (size_t l=L; l>0u; --l, Y-=2u*K, X+=2u*K) { *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
    }

    return 0;
}


int flip_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in flip_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        Y += 2u*(L-1u);
        for (size_t l=L; l>0u; --l, ++X, Y-=2) { *Y = *X; *(Y+1) = *++X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            Y += 2u*(L-1u);
            for (size_t v=V; v>0u; --v, Y+=4u*L)
            {
                for (size_t l=L; l>0u; --l, ++X, Y-=2) { *Y = *X; *(Y+1) = *++X; }
            }
        }
        else if (G==1u)
        {
            Y += 2u*V*(L-1u);
            for (size_t l=L; l>0u; --l, Y-=4u*V)
            {
                for (size_t v=0u; v<2u*V; ++v, ++X, ++Y) { *Y = *X; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y+=2u*K+2u)
                {
                    Y += 2u*K*(L-1u);
                    for (size_t l=L; l>0u; --l, Y-=2u*K, X+=2u*K) { *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
    }

    return 0;
}


int flip_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in flip_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float x1;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L/2u; ++l) { x1 = *X; *X = *(X+L-2u*l-1u); ++X; *(X+L-2u*l-2u) = x1; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=L-L/2u)
            {
                for (size_t l=0u; l<L/2u; ++l) { x1 = *X; *X = *(X+L-2u*l-1u); ++X; *(X+L-2u*l-2u) = x1; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*(L/2u)-1u)
                {
                    for (size_t l=0u; l<L/2u; ++l)
                    {
                        x1 = *X; *X = *(X+K*(L-2u*l-1u));
                        X += K; *(X+K*(L-2u*l-2u)) = x1;
                    }
                }
            }
        }
    }

    return 0;
}


int flip_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in flip_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double x1;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L/2u; ++l) { x1 = *X; *X = *(X+L-2u*l-1u); ++X; *(X+L-2u*l-2u) = x1; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=L-L/2u)
            {
                for (size_t l=0u; l<L/2u; ++l) { x1 = *X; *X = *(X+L-2u*l-1u); ++X; *(X+L-2u*l-2u) = x1; }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*(L/2u)-1u)
                {
                    for (size_t l=0u; l<L/2u; ++l)
                    {
                        x1 = *X; *X = *(X+K*(L-2u*l-1u));
                        X += K; *(X+K*(L-2u*l-2u)) = x1;
                    }
                }
            }
        }
    }

    return 0;
}


int flip_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in flip_inplace_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float x1r, x1i;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L/2u; ++l)
        {
            x1r = *X; x1i = *(X+1);
            *X = *(X+2u*(L-2u*l)-2u); ++X;
            *X = *(X+2u*(L-2u*l)-2u); ++X;
            *(X+2u*(L-2u*l)-4u) = x1r;
            *(X+2u*(L-2u*l)-3u) = x1i;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=2u*(L-L/2u))
            {
                for (size_t l=0u; l<L/2u; ++l)
                {
                    x1r = *X; x1i = *(X+1);
                    *X = *(X+2u*(L-2u*l)-2u); ++X;
                    *X = *(X+2u*(L-2u*l)-2u); ++X;
                    *(X+2u*(L-2u*l)-4u) = x1r;
                    *(X+2u*(L-2u*l)-3u) = x1i;
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*(L/2u)-2u)
                {
                    for (size_t l=0u; l<L/2u; ++l)
                    {
                        x1r = *X; x1i = *(X+1);
                        *X = *(X+2u*K*(L-2u*l-1u)); ++X;
                        *X = *(X+2u*K*(L-2u*l-1u)); X += K-1u;
                        *(X+2u*K*(L-2u*l-2u)) = x1r;
                        *(X+2u*K*(L-2u*l-2u)+1u) = x1i;
                    }
                }
            }
        }
    }

    return 0;
}


int flip_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in flip_inplace_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double x1r, x1i;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<L/2u; ++l)
        {
            x1r = *X; x1i = *(X+1);
            *X = *(X+2u*(L-2u*l)-2u); ++X;
            *X = *(X+2u*(L-2u*l)-2u); ++X;
            *(X+2u*(L-2u*l)-4u) = x1r;
            *(X+2u*(L-2u*l)-3u) = x1i;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, X+=2u*(L-L/2u))
            {
                for (size_t l=0u; l<L/2u; ++l)
                {
                    x1r = *X; x1i = *(X+1);
                    *X = *(X+2u*(L-2u*l)-2u); ++X;
                    *X = *(X+2u*(L-2u*l)-2u); ++X;
                    *(X+2u*(L-2u*l)-4u) = x1r;
                    *(X+2u*(L-2u*l)-3u) = x1i;
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*(L/2u)-2u)
                {
                    for (size_t l=0u; l<L/2u; ++l)
                    {
                        x1r = *X; x1i = *(X+1);
                        *X = *(X+2u*K*(L-2u*l-1u)); ++X;
                        *X = *(X+2u*K*(L-2u*l-1u)); X += K-1u;
                        *(X+2u*K*(L-2u*l-2u)) = x1r;
                        *(X+2u*K*(L-2u*l-2u)+1u) = x1i;
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
