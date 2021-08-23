//Converts radians to degrees for each element of X.
//This has in-place and not-in-place versions.

//The direct for loop was defintely faster that cblas_?scal, even for N=1e6

#include <stdio.h>
//#include <math.h>

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

    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * r2d; }

    return 0;
}


int rad2deg_d (double *Y, const double *X, const size_t N)
{
    const double r2d = 180.0 * M_1_PI;

    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * r2d; }
    
    return 0;
}



int rad2deg_inplace_s (float *X, const size_t N)
{
    const float r2d = (float)(180.0*M_1_PI);

    //cblas_sscal((int)N,r2d,X,1);
    for (size_t n=N; n>0u; --n, ++X) { *X *= r2d; }

    return 0;
}


int rad2deg_inplace_d (double *X, const size_t N)
{
    const double r2d = 180.0 * M_1_PI;

    for (size_t n=N; n>0u; --n, ++X) { *X *= r2d; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
