//Vec2scalar (reduction) operation.
//Gets minimum of values for each vector in X along dim.
//For complex case, finds min absolute value and outputs the complex number.
//For complex case, this is the proper absolute value; see amin for the other one.

//The in-place version was always slower, since requires rewind to start of X,
//so those are no longer included for any of the Vec2scalar functions.
//Instead, in-place will mean for Vec2scalar that X is allowed to be modified.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int min_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int min_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int min_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int min_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int min_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in min_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        mn = *X++;
        for (size_t l=1u; l<L; ++l, ++X) { if (*X<mn) { mn = *X; } }
        *Y = mn;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mn = *X++;
                for (size_t l=1u; l<L; ++l, ++X) { if (*X<mn) { mn = *X; } }
                *Y = mn;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    mn = *X; X += K;
                    for (size_t l=1u; l<L; ++l, X+=K) { if (*X<mn) { mn = *X; } }
                    *Y = mn;
                }
            }
        }
    }

    return 0;
}


int min_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in min_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        mn = *X++;
        for (size_t l=1u; l<L; ++l, ++X) { if (*X<mn) { mn = *X; } }
        *Y = mn;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, ++Y)
            {
                mn = *X++;
                for (size_t l=1u; l<L; ++l, ++X) { if (*X<mn) { mn = *X; } }
                *Y = mn;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1u, ++Y)
                {
                    mn = *X; X += K;
                    for (size_t l=1u; l<L; ++l, X+=K) { if (*X<mn) { mn = *X; } }
                    *Y = mn;
                }
            }
        }
    }

    return 0;
}


int min_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in min_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float xx, mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        mn = *X**X + *(X+1)**(X+1);
        *Y = *X; *(Y+1) = *(X+1); X += 2;
        for (size_t l=1u; l<L; ++l, X+=2)
        {
            xx = *X**X + *(X+1)**(X+1);
            if (xx<mn) { mn = xx; *Y = *X; *(Y+1) = *(X+1); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, Y+=2)
            {
                mn = *X**X + *(X+1)**(X+1);
                *Y = *X; *(Y+1) = *(X+1); X += 2;
                for (size_t l=1u; l<L; ++l, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx<mn) { mn = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, Y+=2)
                {
                    mn = *X**X + *(X+1)**(X+1);
                    *Y = *X; *(Y+1) = *(X+1); X += 2u*K;
                    for (size_t l=1u; l<L; ++l, X+=2u*K)
                    {
                        xx = *X**X + *(X+1)**(X+1);
                        if (xx<mn) { mn = xx; *Y = *X; *(Y+1) = *(X+1); }
                    }
                }
            }
        }
    }

    return 0;
}


int min_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in min_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double xx, mn;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        mn = *X**X + *(X+1)**(X+1);
        *Y = *X; *(Y+1) = *(X+1); X += 2;
        for (size_t l=1u; l<L; ++l, X+=2)
        {
            xx = *X**X + *(X+1)**(X+1);
            if (xx<mn) { mn = xx; *Y = *X; *(Y+1) = *(X+1); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, Y+=2)
            {
                mn = *X**X + *(X+1)**(X+1);
                *Y = *X; *(Y+1) = *(X+1); X += 2;
                for (size_t l=1u; l<L; ++l, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx<mn) { mn = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*L-2u, Y+=2)
                {
                    mn = *X**X + *(X+1)**(X+1);
                    *Y = *X; *(Y+1) = *(X+1); X += 2u*K;
                    for (size_t l=1u; l<L; ++l, X+=2u*K)
                    {
                        xx = *X**X + *(X+1)**(X+1);
                        if (xx<mn) { mn = xx; *Y = *X; *(Y+1) = *(X+1); }
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
