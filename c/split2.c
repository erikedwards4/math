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
    if (dim>3) { fprintf(stderr,"error in split2_s: dim must be in [0 3]\n"); return 1; }
    if (dim==0 && R%2) { fprintf(stderr,"error in split2_s: num rows X must be even for dim=0"); return 1; }
    if (dim==1 && C%2) { fprintf(stderr,"error in split2_s: num cols X must be even for dim=1"); return 1; }
    if (dim==2 && S%2) { fprintf(stderr,"error in split2_s: num slices X must be even for dim=2"); return 1; }
    if (dim==3 && H%2) { fprintf(stderr,"error in split2_s: num hyperslices X must be even for dim=3"); return 1; }

    const size_t B = (iscolmajor) ? ((dim==0) ? R/2 : (dim==1) ? R*C/2 : (dim==2) ? R*C*S/2 : R*C*S*H/2) : ((dim==0) ? H*S*C*R/2 : (dim==1) ? H*S*C/2 : (dim==2) ? H*S/2 : H/2);
    const size_t G = R*C*S*H/(2*B);

    for (size_t g=0; g<G; ++g)
    {
        for (size_t b=0; b<B; ++b, ++X, ++Y1) { *Y1 = *X; }
        for (size_t b=0; b<B; ++b, ++X, ++Y2) { *Y2 = *X; }
        //cblas_scopy((int)B,X,1,Y1,1); X += B; Y1 += B; //slightly faster for some test with very large N
        //cblas_scopy((int)B,X,1,Y2,1); X += B; Y2 += B;
    }

    return 0;
}


int split2_d (double *Y1, double *Y2, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in split2_d: dim must be in [0 3]\n"); return 1; }
    if (dim==0 && R%2) { fprintf(stderr,"error in split2_d: num rows X must be even for dim=0"); return 1; }
    if (dim==1 && C%2) { fprintf(stderr,"error in split2_d: num cols X must be even for dim=1"); return 1; }
    if (dim==2 && S%2) { fprintf(stderr,"error in split2_d: num slices X must be even for dim=2"); return 1; }
    if (dim==3 && H%2) { fprintf(stderr,"error in split2_d: num hyperslices X must be even for dim=3"); return 1; }

    const size_t B = (iscolmajor) ? ((dim==0) ? R/2 : (dim==1) ? R*C/2 : (dim==2) ? R*C*S/2 : R*C*S*H/2) : ((dim==0) ? H*S*C*R/2 : (dim==1) ? H*S*C/2 : (dim==2) ? H*S/2 : H/2);
    const size_t G = R*C*S*H/(2*B);

    for (size_t g=0; g<G; ++g)
    {
        for (size_t b=0; b<B; ++b, ++X, ++Y1) { *Y1 = *X; }
        for (size_t b=0; b<B; ++b, ++X, ++Y2) { *Y2 = *X; }
    }

    return 0;
}


int split2_c (float *Y1, float *Y2, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in split2_c: dim must be in [0 3]\n"); return 1; }
    if (dim==0 && R%2) { fprintf(stderr,"error in split2_c: num rows X must be even for dim=0"); return 1; }
    if (dim==1 && C%2) { fprintf(stderr,"error in split2_c: num cols X must be even for dim=1"); return 1; }
    if (dim==2 && S%2) { fprintf(stderr,"error in split2_c: num slices X must be even for dim=2"); return 1; }
    if (dim==3 && H%2) { fprintf(stderr,"error in split2_c: num hyperslices X must be even for dim=3"); return 1; }

    const size_t B = (iscolmajor) ? ((dim==0) ? R/2 : (dim==1) ? R*C/2 : (dim==2) ? R*C*S/2 : R*C*S*H/2) : ((dim==0) ? H*S*C*R/2 : (dim==1) ? H*S*C/2 : (dim==2) ? H*S/2 : H/2);
    const size_t G = R*C*S*H/(2*B);

    for (size_t g=0; g<G; ++g)
    {
        for (size_t b=0; b<B; ++b, ++X, ++Y1) { *Y1 = *X; *++Y1 = *++X; }
        for (size_t b=0; b<B; ++b, ++X, ++Y2) { *Y2 = *X; *++Y2 = *++X; }
    }

    return 0;
}


int split2_z (double *Y1, double *Y2, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in split2_z: dim must be in [0 3]\n"); return 1; }
    if (dim==0 && R%2) { fprintf(stderr,"error in split2_z: num rows X must be even for dim=0"); return 1; }
    if (dim==1 && C%2) { fprintf(stderr,"error in split2_z: num cols X must be even for dim=1"); return 1; }
    if (dim==2 && S%2) { fprintf(stderr,"error in split2_z: num slices X must be even for dim=2"); return 1; }
    if (dim==3 && H%2) { fprintf(stderr,"error in split2_z: num hyperslices X must be even for dim=3"); return 1; }

    const size_t B = (iscolmajor) ? ((dim==0) ? R/2 : (dim==1) ? R*C/2 : (dim==2) ? R*C*S/2 : R*C*S*H/2) : ((dim==0) ? H*S*C*R/2 : (dim==1) ? H*S*C/2 : (dim==2) ? H*S/2 : H/2);
    const size_t G = R*C*S*H/(2*B);

    for (size_t g=0; g<G; ++g)
    {
        for (size_t b=0; b<B; ++b, ++X, ++Y1) { *Y1 = *X; *++Y1 = *++X; }
        for (size_t b=0; b<B; ++b, ++X, ++Y2) { *Y2 = *X; *++Y2 = *++X; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
