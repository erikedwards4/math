//Vec2scalar (reduction) operation.
//Gets product of elements for each vector in X along dim.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int prod_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int prod_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int prod_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int prod_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int prod_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in prod_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        float pr = *X;
        for (size_t l=1; l<L; ++l, ++X) { pr *= *X; }
        *Y = pr;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            float pr;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                pr = *X++;
                for (size_t l=1; l<L; ++l, ++X) { pr *= *X; }
                *Y = pr;
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y *= *X; }
            }
        }
        else
        {
            float pr;
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    pr = *X; X += K;
                    for (size_t l=1; l<L; ++l, X+=K) { pr *= *X; }
                    *Y = pr;
                }
            }
        }
    }

    return 0;
}


int prod_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in prod_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        double pr = *X;
        for (size_t l=1; l<L; ++l, ++X) { pr *= *X; }
        *Y = pr;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            double pr;
            for (size_t v=0; v<V; ++v, ++Y)
            {
                pr = *X++;
                for (size_t l=1; l<L; ++l, ++X) { pr *= *X; }
                *Y = pr;
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y *= *X; }
            }
        }
        else
        {
            double pr;
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    pr = *X; X += K;
                    for (size_t l=1; l<L; ++l, X+=K) { pr *= *X; }
                    *Y = pr;
                }
            }
        }
    }

    return 0;
}


int prod_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in prod_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float yr, yi;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y = *X++; *(Y+1) = *X++;
        for (size_t l=1; l<L; ++l, X+=2)
        {
            yr = *X**Y - *(X+1)**(Y+1);
            yi = *X**(Y+1) + *(X+1)**Y;
            *Y = yr; *(Y+1) = yi;
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
                *Y = *X++; *(Y+1) = *X++;
                for (size_t l=1; l<L; ++l, X+=2)
                {
                    yr = *X**Y - *(X+1)**(Y+1);
                    yi = *X**(Y+1) + *(X+1)**Y;
                    *Y = yr; *(Y+1) = yi;
                }
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            Y -= 2*V;
            for (size_t l=1; l<L; ++l, Y-=2*V)
            {
                for (size_t v=0; v<V; ++v, X+=2)
                {
                    yr = *X**Y - *(X+1)**(Y+1);
                    yi = *X**(Y+1) + *(X+1)**Y;
                    *Y++ = yr; *Y++ = yi;
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, Y+=2)
                {
                    *Y = *X; *(Y+1) = *(X+1); X += 2*K;
                    for (size_t l=1; l<L; ++l, X+=2*K)
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


int prod_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in prod_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double yr, yi;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *Y = *X++; *(Y+1) = *X++;
        for (size_t l=1; l<L; ++l, X+=2)
        {
            yr = *X**Y - *(X+1)**(Y+1);
            yi = *X**(Y+1) + *(X+1)**Y;
            *Y = yr; *(Y+1) = yi;
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
                *Y = *X++; *(Y+1) = *X++;
                for (size_t l=1; l<L; ++l, X+=2)
                {
                    yr = *X**Y - *(X+1)**(Y+1);
                    yi = *X**(Y+1) + *(X+1)**Y;
                    *Y = yr; *(Y+1) = yi;
                }
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            Y -= 2*V;
            for (size_t l=1; l<L; ++l, Y-=2*V)
            {
                for (size_t v=0; v<V; ++v, X+=2)
                {
                    yr = *X**Y - *(X+1)**(Y+1);
                    yi = *X**(Y+1) + *(X+1)**Y;
                    *Y++ = yr; *Y++ = yi;
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, Y+=2)
                {
                    *Y = *X; *(Y+1) = *(X+1); X += 2*K;
                    for (size_t l=1; l<L; ++l, X+=2*K)
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
