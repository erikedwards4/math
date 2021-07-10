//Vec2scalar (reduction) operation.
//Gets index of minimum of values for each vector in X along dim.
//For complex case, this is the proper absolute value; see iamin for the other one.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int imin_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int imin_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int imin_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int imin_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int imin_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in imin_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        mn = *X++; *Y = 0.0f;
        for (size_t l=1u; l<L; ++l, ++X) { if (*X<mn) { mn = *X; *Y = l; } }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mn = *X++; *Y = 0.0f;
                for (size_t l=1u; l<L; ++l, ++X) { if (*X<mn) { mn = *X; *Y = l; } }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
                {
                    mn = *X; *Y = 0.0f; X += K;
                    for (size_t l=1u; l<L; ++l, X+=K) { if (*X<mn) { mn = *X; *Y = l; } }
                }
            }
        }
    }

    return 0;
}


int imin_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in imin_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        mn = *X++; *Y = 0.0;
        for (size_t l=1u; l<L; ++l, ++X) { if (*X<mn) { mn = *X; *Y = l; } }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mn = *X++; *Y = 0.0;
                for (size_t l=1u; l<L; ++l, ++X) { if (*X<mn) { mn = *X; *Y = l; } }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
                {
                    mn = *X; *Y = 0.0; X += K;
                    for (size_t l=1u; l<L; ++l, X+=K) { if (*X<mn) { mn = *X; *Y = l; } }
                }
            }
        }
    }

    return 0;
}


int imin_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in imin_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float xx, mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    }
    else if (L==N)
    {
        mn = *X**X + *(X+1)**(X+1);
        *Y = 0.0f; X += 2;
        for (size_t l=1u; l<L; ++l, X+=2)
        {
            xx = *X**X + *(X+1)**(X+1);
            if (xx<mn) { mn = xx; *Y = l; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mn = *X**X + *(X+1)**(X+1);
                *Y = 0.0f; X += 2;
                for (size_t l=1u; l<L; ++l, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx<mn) { mn = xx; *Y = l; }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2, ++Y)
                {
                    mn = *X**X + *(X+1)**(X+1);
                    *Y = 0.0f; X += 2*K;
                    for (size_t l=1u; l<L; ++l, X+=2*K)
                    {
                        xx = *X**X + *(X+1)**(X+1);
                        if (xx<mn) { mn = xx; *Y = l; }
                    }
                }
            }
        }
    }

    return 0;
}


int imin_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in imin_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double xx, mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    }
    else if (L==N)
    {
        mn = *X**X + *(X+1)**(X+1);
        *Y = 0.0; X += 2;
        for (size_t l=1u; l<L; ++l, X+=2)
        {
            xx = *X**X + *(X+1)**(X+1);
            if (xx<mn) { mn = xx; *Y = l; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mn = *X**X + *(X+1)**(X+1);
                *Y = 0.0; X += 2;
                for (size_t l=1u; l<L; ++l, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx<mn) { mn = xx; *Y = l; }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2, ++Y)
                {
                    mn = *X**X + *(X+1)**(X+1);
                    *Y = 0.0; X += 2*K;
                    for (size_t l=1u; l<L; ++l, X+=2*K)
                    {
                        xx = *X**X + *(X+1)**(X+1);
                        if (xx<mn) { mn = xx; *Y = l; }
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
