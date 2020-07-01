//Converts radians to degrees for each element of X.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rad2deg_s (float *Y, const float *X, const int N);
int rad2deg_d (double *Y, const double *X, const int N);

int rad2deg_inplace_s (float *X, const int N);
int rad2deg_inplace_d (double *X, const int N);


int rad2deg_s (float *Y, const float *X, const int N)
{
    const float r2d = (float)(180.0*M_1_PI);
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in rad2deg_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = r2d * X[n]; }

    return 0;
}


int rad2deg_d (double *Y, const double *X, const int N)
{
    const double r2d = 180.0 * M_1_PI;
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in rad2deg_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = r2d * X[n]; }
    
    return 0;
}



int rad2deg_inplace_s (float *X, const int N)
{
    const float r2d = (float)(180.0*M_1_PI);

    //Checks
    if (N<0) { fprintf(stderr,"error in rad2deg_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    cblas_sscal(N,r2d,X,1);

    return 0;
}


int rad2deg_inplace_d (double *X, const int N)
{
    const double r2d = 180.0 * M_1_PI;

    //Checks
    if (N<0) { fprintf(stderr,"error in rad2deg_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    cblas_dscal(N,r2d,X,1);
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
