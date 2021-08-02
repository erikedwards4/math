//Splits input X into 3 outputs Y1, Y2, Y3

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int split3_s (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int split3_d (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int split3_c (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int split3_z (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int split3_s (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in split3_s: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R%3u) { fprintf(stderr,"error in split3_s: num rows X must be a multiple of 3 for dim=0"); return 1; }
    if (dim==1u && C%3u) { fprintf(stderr,"error in split3_s: num cols X must be a multiple of 3 for dim=1"); return 1; }
    if (dim==2u && S%3u) { fprintf(stderr,"error in split3_s: num slices X must be a multiple of 3 for dim=2"); return 1; }
    if (dim==3u && H%3u) { fprintf(stderr,"error in split3_s: num hyperslices X must be a multiple of 3 for dim=3"); return 1; }

    const size_t B = (iscolmajor) ? ((dim==0u) ? R/3u : (dim==1u) ? R*C/3u : (dim==2u) ? R*C*S/3u : R*C*S*H/3u) : ((dim==0u) ? H*S*C*R/3u : (dim==1u) ? H*S*C/3u : (dim==2u) ? H*S/3u : H/3u);
    const size_t G = R*C*S*H/(3u*B);

    for (size_t g=0u; g<G; ++g)
    {
        for (size_t b=0u; b<B; ++b, ++X, ++Y1) { *Y1 = *X; }
        for (size_t b=0u; b<B; ++b, ++X, ++Y2) { *Y2 = *X; }
        for (size_t b=0u; b<B; ++b, ++X, ++Y3) { *Y3 = *X; }
    }

    return 0;
}


int split3_d (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in split3_d: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R%3u) { fprintf(stderr,"error in split3_d: num rows X must be a multiple of 3 for dim=0"); return 1; }
    if (dim==1u && C%3u) { fprintf(stderr,"error in split3_d: num cols X must be a multiple of 3 for dim=1"); return 1; }
    if (dim==2u && S%3u) { fprintf(stderr,"error in split3_d: num slices X must be a multiple of 3 for dim=2"); return 1; }
    if (dim==3u && H%3u) { fprintf(stderr,"error in split3_d: num hyperslices X must be a multiple of 3 for dim=3"); return 1; }

    const size_t B = (iscolmajor) ? ((dim==0u) ? R/3u : (dim==1u) ? R*C/3u : (dim==2u) ? R*C*S/3u : R*C*S*H/3u) : ((dim==0u) ? H*S*C*R/3u : (dim==1u) ? H*S*C/3u : (dim==2u) ? H*S/3u : H/3u);
    const size_t G = R*C*S*H/(3u*B);

    for (size_t g=0u; g<G; ++g)
    {
        for (size_t b=0u; b<B; ++b, ++X, ++Y1) { *Y1 = *X; }
        for (size_t b=0u; b<B; ++b, ++X, ++Y2) { *Y2 = *X; }
        for (size_t b=0u; b<B; ++b, ++X, ++Y3) { *Y3 = *X; }
    }

    return 0;
}


int split3_c (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in split3_c: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R%3u) { fprintf(stderr,"error in split3_c: num rows X must be a multiple of 3 for dim=0"); return 1; }
    if (dim==1u && C%3u) { fprintf(stderr,"error in split3_c: num cols X must be a multiple of 3 for dim=1"); return 1; }
    if (dim==2u && S%3u) { fprintf(stderr,"error in split3_c: num slices X must be a multiple of 3 for dim=2"); return 1; }
    if (dim==3u && H%3u) { fprintf(stderr,"error in split3_c: num hyperslices X must be a multiple of 3 for dim=3"); return 1; }

    const size_t B = (iscolmajor) ? ((dim==0u) ? R/3u : (dim==1u) ? R*C/3u : (dim==2u) ? R*C*S/3u : R*C*S*H/3u) : ((dim==0u) ? H*S*C*R/3u : (dim==1u) ? H*S*C/3u : (dim==2u) ? H*S/3u : H/3u);
    const size_t G = R*C*S*H/(3u*B);

    for (size_t g=0u; g<G; ++g)
    {
        for (size_t b=0u; b<B; ++b, ++X, ++Y1) { *Y1 = *X; *++Y1 = *++X; }
        for (size_t b=0u; b<B; ++b, ++X, ++Y2) { *Y2 = *X; *++Y2 = *++X; }
        for (size_t b=0u; b<B; ++b, ++X, ++Y3) { *Y3 = *X; *++Y3 = *++X; }
    }

    return 0;
}


int split3_z (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in split3_z: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R%3u) { fprintf(stderr,"error in split3_z: num rows X must be a multiple of 3 for dim=0"); return 1; }
    if (dim==1u && C%3u) { fprintf(stderr,"error in split3_z: num cols X must be a multiple of 3 for dim=1"); return 1; }
    if (dim==2u && S%3u) { fprintf(stderr,"error in split3_z: num slices X must be a multiple of 3 for dim=2"); return 1; }
    if (dim==3u && H%3u) { fprintf(stderr,"error in split3_z: num hyperslices X must be a multiple of 3 for dim=3"); return 1; }

    const size_t B = (iscolmajor) ? ((dim==0u) ? R/3u : (dim==1u) ? R*C/3u : (dim==2u) ? R*C*S/3u : R*C*S*H/3u) : ((dim==0u) ? H*S*C*R/3u : (dim==1u) ? H*S*C/3u : (dim==2u) ? H*S/3u : H/3u);
    const size_t G = R*C*S*H/(3u*B);

    for (size_t g=0u; g<G; ++g)
    {
        for (size_t b=0u; b<B; ++b, ++X, ++Y1) { *Y1 = *X; *++Y1 = *++X; }
        for (size_t b=0u; b<B; ++b, ++X, ++Y2) { *Y2 = *X; *++Y2 = *++X; }
        for (size_t b=0u; b<B; ++b, ++X, ++Y3) { *Y3 = *X; *++Y3 = *++X; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
