//Splits input X into 2 outputs Y1, Y2

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int split2_s (float *Y1, float *Y2, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int split2_d (double *Y1, double *Y2, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int split2_c (float *Y1, float *Y2, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int split2_z (double *Y1, double *Y2, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int split2_s (float *Y1, float *Y2, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in split2_s: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R%2u) { fprintf(stderr,"error in split2_s: num rows X must be even for dim=0"); return 1; }
    if (dim==1u && C%2u) { fprintf(stderr,"error in split2_s: num cols X must be even for dim=1"); return 1; }
    if (dim==2u && S%2u) { fprintf(stderr,"error in split2_s: num slices X must be even for dim=2"); return 1; }
    if (dim==3u && H%2u) { fprintf(stderr,"error in split2_s: num hyperslices X must be even for dim=3"); return 1; }

    const size_t B = (iscolmajor) ? ((dim==0u) ? R/2u : (dim==1u) ? R*C/2u : (dim==2u) ? R*C*S/2u : R*C*S*H/2u) : ((dim==0u) ? H*S*C*R/2u : (dim==1u) ? H*S*C/2u : (dim==2u) ? H*S/2u : H/2u);
    const size_t G = R*C*S*H/(2u*B);

    for (size_t g=0u; g<G; ++g)
    {
        for (size_t b=0u; b<B; ++b, ++X, ++Y1) { *Y1 = *X; }
        for (size_t b=0u; b<B; ++b, ++X, ++Y2) { *Y2 = *X; }
        //cblas_scopy((int)B,X,1,Y1,1); X += B; Y1 += B; //slightly faster for some test with very large N
        //cblas_scopy((int)B,X,1,Y2,1); X += B; Y2 += B;
    }

    return 0;
}


int split2_d (double *Y1, double *Y2, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in split2_d: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R%2u) { fprintf(stderr,"error in split2_d: num rows X must be even for dim=0"); return 1; }
    if (dim==1u && C%2u) { fprintf(stderr,"error in split2_d: num cols X must be even for dim=1"); return 1; }
    if (dim==2u && S%2u) { fprintf(stderr,"error in split2_d: num slices X must be even for dim=2"); return 1; }
    if (dim==3u && H%2u) { fprintf(stderr,"error in split2_d: num hyperslices X must be even for dim=3"); return 1; }

    const size_t B = (iscolmajor) ? ((dim==0u) ? R/2u : (dim==1u) ? R*C/2u : (dim==2u) ? R*C*S/2u : R*C*S*H/2u) : ((dim==0u) ? H*S*C*R/2u : (dim==1u) ? H*S*C/2u : (dim==2u) ? H*S/2u : H/2u);
    const size_t G = R*C*S*H/(2u*B);

    for (size_t g=0u; g<G; ++g)
    {
        for (size_t b=0u; b<B; ++b, ++X, ++Y1) { *Y1 = *X; }
        for (size_t b=0u; b<B; ++b, ++X, ++Y2) { *Y2 = *X; }
    }

    return 0;
}


int split2_c (float *Y1, float *Y2, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in split2_c: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R%2u) { fprintf(stderr,"error in split2_c: num rows X must be even for dim=0"); return 1; }
    if (dim==1u && C%2u) { fprintf(stderr,"error in split2_c: num cols X must be even for dim=1"); return 1; }
    if (dim==2u && S%2u) { fprintf(stderr,"error in split2_c: num slices X must be even for dim=2"); return 1; }
    if (dim==3u && H%2u) { fprintf(stderr,"error in split2_c: num hyperslices X must be even for dim=3"); return 1; }

    const size_t B = (iscolmajor) ? ((dim==0u) ? R/2u : (dim==1u) ? R*C/2u : (dim==2u) ? R*C*S/2u : R*C*S*H/2u) : ((dim==0u) ? H*S*C*R/2u : (dim==1u) ? H*S*C/2u : (dim==2u) ? H*S/2u : H/2u);
    const size_t G = R*C*S*H/(2u*B);

    for (size_t g=0u; g<G; ++g)
    {
        for (size_t b=0u; b<B; ++b, ++X, ++Y1) { *Y1 = *X; *++Y1 = *++X; }
        for (size_t b=0u; b<B; ++b, ++X, ++Y2) { *Y2 = *X; *++Y2 = *++X; }
    }

    return 0;
}


int split2_z (double *Y1, double *Y2, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in split2_z: dim must be in [0 3]\n"); return 1; }
    if (dim==0u && R%2u) { fprintf(stderr,"error in split2_z: num rows X must be even for dim=0"); return 1; }
    if (dim==1u && C%2u) { fprintf(stderr,"error in split2_z: num cols X must be even for dim=1"); return 1; }
    if (dim==2u && S%2u) { fprintf(stderr,"error in split2_z: num slices X must be even for dim=2"); return 1; }
    if (dim==3u && H%2u) { fprintf(stderr,"error in split2_z: num hyperslices X must be even for dim=3"); return 1; }

    const size_t B = (iscolmajor) ? ((dim==0u) ? R/2u : (dim==1u) ? R*C/2u : (dim==2u) ? R*C*S/2u : R*C*S*H/2u) : ((dim==0u) ? H*S*C*R/2u : (dim==1u) ? H*S*C/2u : (dim==2u) ? H*S/2u : H/2u);
    const size_t G = R*C*S*H/(2u*B);

    for (size_t g=0u; g<G; ++g)
    {
        for (size_t b=0u; b<B; ++b, ++X, ++Y1) { *Y1 = *X; *++Y1 = *++X; }
        for (size_t b=0u; b<B; ++b, ++X, ++Y2) { *Y2 = *X; *++Y2 = *++X; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
