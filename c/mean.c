//Vec2scalar (reduction) operation.
//Gets mean for each vector in X along dim.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

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
    const float ni = 1.0f / L;

    if (N==0) {}
    else if (L==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        if (L<7000)
        {
            *Y = 0.0f;
            for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
            *Y *= ni;
        }
        else
        {
            const float o = 1.0f;
            *Y = cblas_sdot((int)L,X,1,&o,0) * ni;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            if (L<80)
            {
                for (size_t v=0; v<V; ++v, ++Y)
                {
                    *Y = 0.0f;
                    for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
                    *Y *= ni;
                }
            }
            else
            {
                const float o = 1.0f;
                for (size_t v=0; v<V; ++v, X+=L, ++Y)
                {
                    *Y = cblas_sdot((int)L,X,1,&o,0) * ni;
                }
            }
        }
        else if (G==1)
        {
            const float z = 0.0f;
            cblas_scopy((int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += *X; }
            }
            cblas_sscal((int)V,ni,Y,1);
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=L*K-1, ++Y)
                {
                    *Y = 0.0f;
                    for (size_t l=0; l<L; ++l, X+=K) { *Y += *X; }
                    *Y *= ni;
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
    const double ni = 1.0 / L;

    if (N==0) {}
    else if (L==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        if (L<7000)
        {
            *Y = 0.0;
            for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
            *Y *= ni;
        }
        else
        {
            const double o = 1.0;
            *Y = cblas_ddot((int)L,X,1,&o,0) * ni;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            if (L<80)
            {
                for (size_t v=0; v<V; ++v, ++Y)
                {
                    *Y = 0.0;
                    for (size_t l=0; l<L; ++l, ++X) { *Y += *X; }
                    *Y *= ni;
                }
            }
            else
            {
                const double o = 1.0;
                for (size_t v=0; v<V; ++v, X+=L, ++Y)
                {
                    *Y = cblas_ddot((int)L,X,1,&o,0) * ni;
                }
            }
        }
        else if (G==1)
        {
            const double z = 0.0;
            cblas_dcopy((int)V,&z,0,Y,1);
            for (size_t l=0; l<L; ++l, Y-=V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y += *X; }
            }
            cblas_dscal((int)V,ni,Y,1);
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=L*K-1, ++Y)
                {
                    *Y = 0.0;
                    for (size_t l=0; l<L; ++l, X+=K) { *Y += *X; }
                    *Y *= ni;
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
    const float ni = 1.0f / L;
    float yr, yi;

    if (N==0) {}
    else if (L==1) { cblas_ccopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        yr = *X++; yi = *X++;
        for (size_t l=1; l<L; ++l, ++X) { yr += *X++; yi += *X; }
        *Y++ = yr * ni; *Y = yi *ni;
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
                yr = *X++; yi = *X++;
                for (size_t l=1; l<L; ++l, ++X) { yr += *X++; yi += *X; }
                *Y++ = yr * ni; *Y = yi * ni;
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y++ = *X++; *Y = *X; }
            Y -= 2*V;
            for (size_t l=1; l<L; ++l, Y-=2*V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y++ += *X++; *Y += *X; }
            }
            cblas_csscal((int)V,ni,Y,1);
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*L*K-2, ++Y)
                {
                    yr = *X++; yi = *X; X += 2*K-1;
                    for (size_t l=1; l<L; ++l, X+=2*K-1) { yr += *X++; yi += *X; }
                    *Y++ = yr * ni; *Y = yi * ni;
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
    const double ni = 1.0 / L;
    double yr, yi;

    if (N==0) {}
    else if (L==1) { cblas_zcopy((int)N,X,1,Y,1); }
    else if (L==N)
    {
        yr = *X++; yi = *X++;
        for (size_t l=1; l<L; ++l, ++X) { yr += *X++; yi += *X; }
        *Y++ = yr * ni; *Y = yi *ni;
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
                yr = *X++; yi = *X++;
                for (size_t l=1; l<L; ++l, ++X) { yr += *X++; yi += *X; }
                *Y++ = yr * ni; *Y = yi * ni;
            }
        }
        else if (G==1)
        {
            for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y++ = *X++; *Y = *X; }
            Y -= 2*V;
            for (size_t l=1; l<L; ++l, Y-=2*V)
            {
                for (size_t v=0; v<V; ++v, ++X, ++Y) { *Y++ += *X++; *Y += *X; }
            }
            cblas_zdscal((int)V,ni,Y,1);
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*L*K-2, ++Y)
                {
                    yr = *X++; yi = *X; X += 2*K-1;
                    for (size_t l=1; l<L; ++l, X+=2*K-1) { yr += *X++; yi += *X; }
                    *Y++ = yr * ni; *Y = yi * ni;
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
