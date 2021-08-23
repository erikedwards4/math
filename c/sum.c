//Vec2scalar (reduction) operation.
//Gets sum of elements for each vector in X along dim.
//For complex case, real and imag parts calculated separately.

//For 2D case, this is definitely faster using cblas_?gemv.

//No in-place version, since cblas_?gemv doesn't work for that.
//Also, I decided against in-place versions in general for vec2scalar operations,
//since this would require rewinding X to the right element for each result.

//I could not confirm any speed gain for cblas_?dot method; much slower for most tests.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sum_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int sum_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int sum_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int sum_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int sum_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in sum_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (L<100u)
        {
            *Y = 0.0f;
            for (size_t l=L; l>0u; --l, ++X) { *Y += *X; }
        }
        else
        {
            float sm = 0.0f;
            for (size_t l=L; l>0u; --l, ++X) { sm += *X; }
            *Y = sm;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (L<10u)
            {
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    *Y = 0.0f;
                    for (size_t l=L; l>0u; --l, ++X) { *Y += *X; }
                }
            }
            else
            {
                float sm;
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=L; l>0u; --l, ++X) { sm += *X; }
                    *Y = sm;
                }
            }
        }
        else if (G==1u)
        {
            //float sm[V]; //this is only C99 allowed, and not faster!
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += *X; }
            }
        }
        else
        {
            float sm;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=K) { *Y += *X; }
                    *Y = sm;
                }
            }
        }
    }

    return 0;
}


int sum_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in sum_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        if (L<100u)
        {
            *Y = 0.0;
            for (size_t l=L; l>0u; --l, ++X) { *Y += *X; }
        }
        else
        {
            double sm = 0.0;
            for (size_t l=L; l>0u; --l, ++X) { sm += *X; }
            *Y = sm;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (L<10u)
            {
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    *Y = 0.0;
                    for (size_t l=L; l>0u; --l, ++X) { *Y += *X; }
                }
            }
            else
            {
                double sm;
                for (size_t v=V; v>0u; --v, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=L; l>0u; --l, ++X) { sm += *X; }
                    *Y = sm;
                }
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; }
            Y -= V;
            for (size_t l=1u; l<L; ++l, Y-=V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += *X; }
            }
        }
        else
        {
            double sm;
            for (size_t g=G; g>0u; --g, X+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-1u, ++Y)
                {
                    sm = 0.0;
                    for (size_t l=L; l>0u; --l, X+=K) { *Y += *X; }
                    *Y = sm;
                }
            }
        }
    }

    return 0;
}


int sum_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in sum_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        float yr = 0.0f, yi = 0.0f;
        for (size_t l=L; l>0u; --l, ++X) { yr += *X; yi += *++X; }
        *Y = yr; *++Y = yi;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            float yr, yi;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                yr = yi = 0.0f;
                for (size_t l=L; l>0u; --l, ++X) { yr += *X; yi += *++X; }
                *Y = yr; *++Y = yi;
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            Y -= 2u*V;
            for (size_t l=1u; l<L; ++l, Y-=2u*V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += *X; *++Y += *++X; }
            }
        }
        else
        {
            float yr, yi;
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    yr = yi = 0.0f;
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u) { yr += *X; yi += *++X; }
                    *Y = yr; *++Y = yi;
                }
            }
        }
    }

    return 0;
}


int sum_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in sum_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        double yr = 0.0, yi = 0.0;
        for (size_t l=L; l>0u; --l, ++X) { yr += *X; yi += *++X; }
        *Y = yr; *++Y = yi;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            double yr, yi;
            for (size_t v=V; v>0u; --v, ++Y)
            {
                yr = yi = 0.0;
                for (size_t l=L; l>0u; --l, ++X) { yr += *X; yi += *++X; }
                *Y = yr; *++Y = yi;
            }
        }
        else if (G==1u)
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y = *X; *++Y = *++X; }
            Y -= 2u*V;
            for (size_t l=1u; l<L; ++l, Y-=2u*V)
            {
                for (size_t v=V; v>0u; --v, ++X, ++Y) { *Y += *X; *++Y += *++X; }
            }
        }
        else
        {
            double yr, yi;
            for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, ++Y)
                {
                    yr = yi = 0.0;
                    for (size_t l=L; l>0u; --l, X+=2u*K-1u) { yr += *X; yi += *++X; }
                    *Y = yr; *++Y = yi;
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
