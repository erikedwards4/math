//Vec2scalar (reduction) operation.
//Gets index of minimum of values for each vector in X along dim.
//For complex case, this is the proper absolute value; see iamin for the other one.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int imin_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int imin_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int imin_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int imin_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int imin_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in imin_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        mn = *X++; *Y = 0.0f;
        for (size_t l=1u; l<L; ++l, ++X) { if (*X<mn) { mn = *X; *Y = (float)l; } }
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
                mn = *X++; *Y = 0.0f;
                for (size_t l=1u; l<L; ++l, ++X) { if (*X<mn) { mn = *X; *Y = (float)l; } }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    mn = *X; *Y = 0.0f; X += K;
                    for (size_t l=1u; l<L; ++l, X+=K) { if (*X<mn) { mn = *X; *Y = (float)l; } }
                }
            }
        }
    }

    return 0;
}


int imin_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in imin_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        mn = *X++; *Y = 0.0;
        for (size_t l=1u; l<L; ++l, ++X) { if (*X<mn) { mn = *X; *Y = (double)l; } }
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
                mn = *X++; *Y = 0.0;
                for (size_t l=1u; l<L; ++l, ++X) { if (*X<mn) { mn = *X; *Y = (double)l; } }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    mn = *X; *Y = 0.0; X += K;
                    for (size_t l=1u; l<L; ++l, X+=K) { if (*X<mn) { mn = *X; *Y = (double)l; } }
                }
            }
        }
    }

    return 0;
}


int imin_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in imin_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float xx, mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        mn = *X**X + *(X+1)**(X+1);
        *Y = 0.0f; X += 2;
        for (size_t l=1u; l<L; ++l, X+=2)
        {
            xx = *X**X + *(X+1)**(X+1);
            if (xx<mn) { mn = xx; *Y = (float)l; }
        }
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
                mn = *X**X + *(X+1)**(X+1);
                *Y = 0.0f; X += 2;
                for (size_t l=1u; l<L; ++l, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx<mn) { mn = xx; *Y = (float)l; }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    mn = *X**X + *(X+1)**(X+1);
                    *Y = 0.0f; X += 2u*K;
                    for (size_t l=1u; l<L; ++l, X+=2u*K)
                    {
                        xx = *X**X + *(X+1)**(X+1);
                        if (xx<mn) { mn = xx; *Y = (float)l; }
                    }
                }
            }
        }
    }

    return 0;
}


int imin_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in imin_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double xx, mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        mn = *X**X + *(X+1)**(X+1);
        *Y = 0.0; X += 2;
        for (size_t l=1u; l<L; ++l, X+=2)
        {
            xx = *X**X + *(X+1)**(X+1);
            if (xx<mn) { mn = xx; *Y = (double)l; }
        }
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
                mn = *X**X + *(X+1)**(X+1);
                *Y = 0.0; X += 2;
                for (size_t l=1u; l<L; ++l, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx<mn) { mn = xx; *Y = (double)l; }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    mn = *X**X + *(X+1)**(X+1);
                    *Y = 0.0; X += 2u*K;
                    for (size_t l=1u; l<L; ++l, X+=2u*K)
                    {
                        xx = *X**X + *(X+1)**(X+1);
                        if (xx<mn) { mn = xx; *Y = (double)l; }
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
