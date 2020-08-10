//Vec2scalar (reduction) operation.
//Gets mean for each vector in X along dim.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int mean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int mean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / L;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (L<100)
        {
            *Y = 0.0f;
            for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
            *Y *= den;
        }
        else
        {
            float sm = 0.0f;
            for (size_t l=0; l<L; ++l, ++X) { sm += *X; }
            *Y = sm * den;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            if (L<10)
            {
                for (size_t v=0; v<V; ++v, ++Y)
                {
                    *Y = 0.0f;
                    for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
                    *Y *= den;
                }
            }
            else
            {
                float sm;
                for (size_t v=0; v<V; ++v, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
                    *Y = sm * den;
                }
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += *X; }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y *= den; }
        }
        else
        {
            float sm;
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=0; l<L; ++l, X+=K) { sm += *X; }
                    *Y = sm * den;
                }
            }
        }
    }

    return 0;
}


int mean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / L;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (L<100)
        {
            *Y = 0.0;
            for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
            *Y *= den;
        }
        else
        {
            double sm = 0.0;
            for (size_t l=0; l<L; ++l, ++X) { sm += *X; }
            *Y = sm * den;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            if (L<10)
            {
                for (size_t v=0; v<V; ++v, ++Y)
                {
                    *Y = 0.0;
                    for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
                    *Y *= den;
                }
            }
            else
            {
                double sm;
                for (size_t v=0; v<V; ++v, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
                    *Y = sm * den;
                }
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += *X; }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y *= den; }
        }
        else
        {
            double sm;
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=0; l<L; ++l, X+=K) { sm += *X; }
                    *Y = sm * den;
                }
            }
        }
    }

    return 0;
}


int mean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float den = 1.0f / L;
    float yr, yi;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        yr = yi = 0.0f;
        for (size_t l=0; l<L; ++l, ++X) { yr += *X; yi += *++X; }
        *Y = yr * den; *++Y = yi * den;
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
                yr = yi = 0.0f;
                for (size_t l=0; l<L; ++l, ++X) { yr += *X; yi += *++X; }
                *Y = yr * den; *++Y = yi * den;
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            Y -= 2*V;
            for (size_t l=1; l<L; ++l, Y-=2*V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += *X; *++Y += *++X; }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y *= den; *++Y *= den; }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, ++Y)
                {
                    yr = yi = 0.0f;
                    for (size_t l=0; l<L; ++l, X+=2*K-1) { yr += *X; yi += *++X; }
                    *Y = yr * den; *++Y = yi * den;
                }
            }
        }
    }

    return 0;
}


int mean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in mean_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double den = 1.0 / L;
    double yr, yi;

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        yr = yi = 0.0;
        for (size_t l=0; l<L; ++l, ++X) { yr += *X; yi += *++X; }
        *Y = yr * den; *++Y = yi * den;
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
                yr = yi = 0.0;
                for (size_t l=0; l<L; ++l, ++X) { yr += *X; yi += *++X; }
                *Y = yr * den; *++Y = yi * den;
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            Y -= 2*V;
            for (size_t l=1; l<L; ++l, Y-=2*V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += *X; *++Y += *++X; }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y *= den; *++Y *= den; }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, ++Y)
                {
                    yr = yi = 0.0;
                    for (size_t l=0; l<L; ++l, X+=2*K-1) { yr += *X; yi += *++X; }
                    *Y = yr * den; *++Y = yi * den;
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
