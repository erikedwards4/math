//Splits input X into 3 outputs Y1, Y2, Y3

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int split3_s (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int split3_d (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int split3_c (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int split3_z (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int split3_s (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in split3_s: dim must be in [0 3]\n"); return 1; }
    if (dim==0 && R%3) { fprintf(stderr,"error in split3_s: num rows X must be a multiple of 3 for dim=0"); return 1; }
    if (dim==1 && C%3) { fprintf(stderr,"error in split3_s: num cols X must be a multiple of 3 for dim=1"); return 1; }
    if (dim==2 && S%3) { fprintf(stderr,"error in split3_s: num slices X must be a multiple of 3 for dim=2"); return 1; }
    if (dim==3 && H%3) { fprintf(stderr,"error in split3_s: num hyperslices X must be a multiple of 3 for dim=3"); return 1; }

    const size_t M = (iscolmajor) ? ((dim==0) ? R/3 : (dim==1) ? R*C/3 : (dim==2) ? R*C*S/3 : R*C*S*H/3) : ((dim==0) ? H*S*C*R/3 : (dim==1) ? H*S*C/3 : (dim==2) ? H*S/3 : H/3);
    const size_t L = R*C*S*H/(3*M);

    for (size_t l=0; l<L; l++)
    {
        cblas_scopy((int)M,X,1,Y1,1); X += M; Y1 += M; 
        cblas_scopy((int)M,X,1,Y2,1); X += M; Y2 += M;
        cblas_scopy((int)M,X,1,Y3,1); X += M; Y3 += M;
    }

    return 0;
}


int split3_d (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in split3_d: dim must be in [0 3]\n"); return 1; }
    if (dim==0 && R%3) { fprintf(stderr,"error in split3_d: num rows X must be a multiple of 3 for dim=0"); return 1; }
    if (dim==1 && C%3) { fprintf(stderr,"error in split3_d: num cols X must be a multiple of 3 for dim=1"); return 1; }
    if (dim==2 && S%3) { fprintf(stderr,"error in split3_d: num slices X must be a multiple of 3 for dim=2"); return 1; }
    if (dim==3 && H%3) { fprintf(stderr,"error in split3_d: num hyperslices X must be a multiple of 3 for dim=3"); return 1; }

    const size_t M = (iscolmajor) ? ((dim==0) ? R/3 : (dim==1) ? R*C/3 : (dim==2) ? R*C*S/3 : R*C*S*H/3) : ((dim==0) ? H*S*C*R/3 : (dim==1) ? H*S*C/3 : (dim==2) ? H*S/3 : H/3);
    const size_t L = R*C*S*H/(3*M);

    for (size_t l=0; l<L; l++)
    {
        cblas_dcopy((int)M,X,1,Y1,1); X += M; Y1 += M; 
        cblas_dcopy((int)M,X,1,Y2,1); X += M; Y2 += M;
        cblas_dcopy((int)M,X,1,Y3,1); X += M; Y3 += M;
    }

    return 0;
}


int split3_c (float *Y1, float *Y2, float *Y3, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in split3_c: dim must be in [0 3]\n"); return 1; }
    if (dim==0 && R%3) { fprintf(stderr,"error in split3_c: num rows X must be a multiple of 3 for dim=0"); return 1; }
    if (dim==1 && C%3) { fprintf(stderr,"error in split3_c: num cols X must be a multiple of 3 for dim=1"); return 1; }
    if (dim==2 && S%3) { fprintf(stderr,"error in split3_c: num slices X must be a multiple of 3 for dim=2"); return 1; }
    if (dim==3 && H%3) { fprintf(stderr,"error in split3_c: num hyperslices X must be a multiple of 3 for dim=3"); return 1; }

    const size_t M = (iscolmajor) ? ((dim==0) ? R/3 : (dim==1) ? R*C/3 : (dim==2) ? R*C*S/3 : R*C*S*H/3) : ((dim==0) ? H*S*C*R/3 : (dim==1) ? H*S*C/3 : (dim==2) ? H*S/3 : H/3);
    const size_t L = R*C*S*H/(3*M);

    for (size_t l=0; l<L; l++)
    {
        cblas_ccopy((int)M,X,1,Y1,1); X += M; Y1 += M; 
        cblas_ccopy((int)M,X,1,Y2,1); X += M; Y2 += M;
        cblas_ccopy((int)M,X,1,Y3,1); X += M; Y3 += M;
    }

    return 0;
}


int split3_z (double *Y1, double *Y2, double *Y3, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in split3_z: dim must be in [0 3]\n"); return 1; }
    if (dim==0 && R%3) { fprintf(stderr,"error in split3_z: num rows X must be a multiple of 3 for dim=0"); return 1; }
    if (dim==1 && C%3) { fprintf(stderr,"error in split3_z: num cols X must be a multiple of 3 for dim=1"); return 1; }
    if (dim==2 && S%3) { fprintf(stderr,"error in split3_z: num slices X must be a multiple of 3 for dim=2"); return 1; }
    if (dim==3 && H%3) { fprintf(stderr,"error in split3_z: num hyperslices X must be a multiple of 3 for dim=3"); return 1; }

    const size_t M = (iscolmajor) ? ((dim==0) ? R/3 : (dim==1) ? R*C/3 : (dim==2) ? R*C*S/3 : R*C*S*H/3) : ((dim==0) ? H*S*C*R/3 : (dim==1) ? H*S*C/3 : (dim==2) ? H*S/3 : H/3);
    const size_t L = R*C*S*H/(3*M);

    for (size_t l=0; l<L; l++)
    {
        cblas_zcopy((int)M,X,1,Y1,1); X += M; Y1 += M; 
        cblas_zcopy((int)M,X,1,Y2,1); X += M; Y2 += M;
        cblas_zcopy((int)M,X,1,Y3,1); X += M; Y3 += M;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
