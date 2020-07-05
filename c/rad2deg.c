//Converts radians to degrees for each element of X.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifndef M_1_PI
    #define M_1_PI 0.318309886183790671538
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rad2deg_s (float *Y, const float *X, const size_t N);
int rad2deg_d (double *Y, const double *X, const size_t N);

int rad2deg_inplace_s (float *X, const size_t N);
int rad2deg_inplace_d (double *X, const size_t N);


int rad2deg_s (float *Y, const float *X, const size_t N)
{
    const float r2d = (float)(180.0*M_1_PI);

    for (size_t n=0; n<N; n++) { Y[n] = r2d * X[n]; }

    return 0;
}


int rad2deg_d (double *Y, const double *X, const size_t N)
{
    const double r2d = 180.0 * M_1_PI;

    for (size_t n=0; n<N; n++) { Y[n] = r2d * X[n]; }
    
    return 0;
}



int rad2deg_inplace_s (float *X, const size_t N)
{
    const float r2d = (float)(180.0*M_1_PI);

    cblas_sscal((int)N,r2d,X,1);

    return 0;
}


int rad2deg_inplace_d (double *X, const size_t N)
{
    const double r2d = 180.0 * M_1_PI;

    cblas_dscal((int)N,r2d,X,1);
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
