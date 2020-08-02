//Vec2scalar (reduction) operation.
//Gets harmonic mean for each vector in X along dim.
//This is the reciprocal of the mean of reciprocals for each vector.

//The same definition is used for complex numbers (as in Octave),
//but there are some publications on this topic, so check into later if needed.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int harmean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int harmean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int harmean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int harmean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int harmean_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in harmean_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y = 0.0f;
        for (size_t l=0; l<L; ++l, ++X) { *Y += 1.0f / *X; }
        *Y = L / *Y;
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
                *Y = 0.0f;
                for (size_t l=0; l<L; ++l, ++X) { *Y += 1.0f / *X; }
                *Y = L / *Y;
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = 1.0f / *X; }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += 1.0f / *X; }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = L / *Y; }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    *Y = 0.0f;
                    for (size_t l=0; l<L; ++l, X+=K) { *Y += 1.0f / *X; }
                    *Y = L / *Y;
                }
            }
        }
    }

    return 0;
}


int harmean_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in harmean_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==0) {}
    else if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        *Y = 0.0;
        for (size_t l=0; l<L; ++l, ++X) { *Y += 1.0 / *X; }
        *Y = L / *Y;
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
                *Y = 0.0;
                for (size_t l=0; l<L; ++l, ++X) { *Y += 1.0 / *X; }
                *Y = L / *Y;
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y = 1.0 / *X; }
            Y -= V;
            for (size_t l=1; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += 1.0 / *X; }
            }
            for (size_t v=0; v<V; ++v, ++Y) { *Y = L / *Y; }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, ++Y)
                {
                    *Y = 0.0;
                    for (size_t l=0; l<L; ++l, X+=K) { *Y += 1.0 / *X; }
                    *Y = L / *Y;
                }
            }
        }
    }

    return 0;
}


int harmean_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in harmean_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float xr, sq, yr, yi;

    if (N==0) {}
    else if (L==1) { cblas_ccopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        yr = yi = 0.0f;
        for (size_t l=0; l<L; ++l, ++X)
        {
            xr = *X++;
            sq = xr*xr + *X**X;
            yr += xr/sq; yi += *X/sq;
        }
        sq = yr*yr + yi*yi;
        *Y++ = yr*L/sq; *Y = yi*L/sq;
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
                for (size_t l=0; l<L; ++l, ++X)
                {
                    xr = *X++;
                    sq = xr*xr + *X**X;
                    yr += xr/sq; yi += *X/sq;
                }
                sq = yr*yr + yi*yi;
                *Y++ = yr*L/sq; *Y = yi*L/sq;
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, ++Y)
                {
                    yr = yi = 0.0f;
                    for (size_t l=0; l<L; ++l, X+=2*K-1)
                    {
                        xr = *X++;
                        sq = xr*xr + *X**X;
                        yr += xr/sq; yi += *X/sq;
                    }
                    sq = yr*yr + yi*yi;
                    *Y++ = yr*L/sq; *Y = yi*L/sq;
                }
            }
        }
    }

    return 0;
}


int harmean_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in harmean_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double xr, sq, yr, yi;

    if (N==0) {}
    else if (L==1) { cblas_zcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        yr = yi = 0.0;
        for (size_t l=0; l<L; ++l, ++X)
        {
            xr = *X++;
            sq = xr*xr + *X**X;
            yr += xr/sq; yi += *X/sq;
        }
        sq = yr*yr + yi*yi;
        *Y++ = yr*L/sq; *Y = yi*L/sq;
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
                for (size_t l=0; l<L; ++l, ++X)
                {
                    xr = *X++;
                    sq = xr*xr + *X**X;
                    yr += xr/sq; yi += *X/sq;
                }
                sq = yr*yr + yi*yi;
                *Y++ = yr*L/sq; *Y = yi*L/sq;
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, ++Y)
                {
                    yr = yi = 0.0;
                    for (size_t l=0; l<L; ++l, X+=2*K-1)
                    {
                        xr = *X++;
                        sq = xr*xr + *X**X;
                        yr += xr/sq; yi += *X/sq;
                    }
                    sq = yr*yr + yi*yi;
                    *Y++ = yr*L/sq; *Y = yi*L/sq;
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
