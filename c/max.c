//Vec2scalar (reduction) operation.
//Gets maximum of values for each vector in X along dim.
//For complex case, finds max absolute value and outputs the complex number.
//For complex case, this is the proper absolute value; see amax for the other one.

#include <stdio.h>
#include <stdlib.h>

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

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y = *X++;
        for (size_t l=1; l<L; ++l, ++X) { if (*X>*Y) { *Y = *X; } }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, ++Y)
            {
                *Y = *X++;
                for (size_t l=1; l<L; ++l, ++X) { if (*X>*Y) { *Y = *X; } }
            }
        }
        else if (G==1)
        {
            for (size_t l=0; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y)
                {
                    if (l==0 || *X>*Y) { *Y = *X; }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    *Y = *X; X += K;
                    for (size_t l=1; l<L; ++l, X+=K) { if (*X>*Y) { *Y = *X; } }
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

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y = *X++;
        for (size_t l=1; l<L; ++l, ++X) { if (*X>*Y) { *Y = *X; } }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v, ++Y)
            {
                *Y = *X++;
                for (size_t l=1; l<L; ++l, ++X) { if (*X>*Y) { *Y = *X; } }
            }
        }
        else if (G==1)
        {
            for (size_t l=0; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y)
                {
                    if (l==0 || *X>*Y) { *Y = *X; }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    *Y = *X; X += K;
                    for (size_t l=1; l<L; ++l, X+=K) { if (*X>*Y) { *Y = *X; } }
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

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        mx = *X**X + *(X+1)**(X+1);
        *Y = *X; *(Y+1) = *(X+1); X += 2;
        for (size_t l=1; l<L; ++l, X+=2)
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
            for (size_t v=0; v<V; ++v, Y+=2)
            {
                mx = *X**X + *(X+1)**(X+1);
                *Y = *X; *(Y+1) = *(X+1); X += 2;
                for (size_t l=1; l<L; ++l, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx>mx) { mx = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
        else if (G==1)
        {
            float *mxs;
            if (!(mxs=(float *)malloc(V*sizeof(float)))) { fprintf(stderr,"error in max_c: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t l=0; l<L; ++l, Y-=2*V, mxs-=V)
            {
                for (size_t v=0; v<V; ++v, X+=2, Y+=2, ++mxs)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (l==0 || xx>*mxs) { *mxs = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
            free(mxs);
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, Y+=2)
                {
                    mx = *X**X + *(X+1)**(X+1);
                    *Y = *X; *(Y+1) = *(X+1); X += 2*K;
                    for (size_t l=1; l<L; ++l, X+=2*K)
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

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        mx = *X**X + *(X+1)**(X+1);
        *Y = *X; *(Y+1) = *(X+1); X += 2;
        for (size_t l=1; l<L; ++l, X+=2)
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
            for (size_t v=0; v<V; ++v, Y+=2)
            {
                mx = *X**X + *(X+1)**(X+1);
                *Y = *X; *(Y+1) = *(X+1); X += 2;
                for (size_t l=1; l<L; ++l, X+=2)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (xx>mx) { mx = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
        else if (G==1)
        {
            double *mxs;
            if (!(mxs=(double *)malloc(V*sizeof(double)))) { fprintf(stderr,"error in max_z: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t l=0; l<L; ++l, Y-=2*V, mxs-=V)
            {
                for (size_t v=0; v<V; ++v, X+=2, Y+=2, ++mxs)
                {
                    xx = *X**X + *(X+1)**(X+1);
                    if (l==0 || xx>*mxs) { *mxs = xx; *Y = *X; *(Y+1) = *(X+1); }
                }
            }
            free(mxs);
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, Y+=2)
                {
                    mx = *X**X + *(X+1)**(X+1);
                    *Y = *X; *(Y+1) = *(X+1); X += 2*K;
                    for (size_t l=1; l<L; ++l, X+=2*K)
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
