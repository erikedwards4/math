//Vec2scalar (reduction) operation.
//Gets maximum of values for each vector in X along dim.
//For complex case, finds max absolute value and outputs the complex number.
//For complex case, this is the proper absolute value; see amax for the other one.

//Computational note: it is much faster to use intermediate mx variable (on stack) than to compare to *Y.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int max_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int max_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int max_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int max_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int max_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in max_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float mx;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        mx = *X++;
        for (size_t l=1u; l<L; ++l, ++X) { if (*X>mx) { mx = *X; } }
        *Y = mx;
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
                mx = *X++;
                for (size_t l=1u; l<L; ++l, ++X) { if (*X>mx) { mx = *X; } }
                *Y = mx;
            }
        }
        // else if (G==1)
        // {
        //     for (size_t l=0u; l<L; ++l, Y-=V)
        //     {
        //         for (size_t v=0u; v<V; ++v, ++X, ++Y)
        //         {
        //             if (l==0 || *X>*Y) { *Y = *X; }
        //         }
        //     }
        // }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
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


int max_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in max_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double mx;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        mx = *X++;
        for (size_t l=1u; l<L; ++l, ++X) { if (*X>mx) { mx = *X; } }
        *Y = mx;
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
                mx = *X++;
                for (size_t l=1u; l<L; ++l, ++X) { if (*X>mx) { mx = *X; } }
                *Y = mx;
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=K*L-1, ++Y)
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


int max_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in max_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float xx, mx;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
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
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, Y+=2)
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
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2, Y+=2)
                {
                    mx = *X**X + *(X+1)**(X+1);
                    *Y = *X; *(Y+1) = *(X+1); X += 2*K;
                    for (size_t l=1u; l<L; ++l, X+=2*K)
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


int max_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in max_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double xx, mx;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=0u; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
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
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0u; v<V; ++v, Y+=2)
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
            for (size_t g=0u; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0u; b<B; ++b, X-=2*K*L-2, Y+=2)
                {
                    mx = *X**X + *(X+1)**(X+1);
                    *Y = *X; *(Y+1) = *(X+1); X += 2*K;
                    for (size_t l=1u; l<L; ++l, X+=2*K)
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
