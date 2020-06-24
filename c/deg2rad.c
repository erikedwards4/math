//Converts degrees to radians for each element of X.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int deg2rad_s (float *Y, const float *X, const int N);
int deg2rad_d (double *Y, const double *X, const int N);

int deg2rad_inplace_s (float *X, const int N);
int deg2rad_inplace_d (double *X, const int N);


int deg2rad_s (float *Y, const float *X, const int N)
{
    const float d2r = (float)(M_PI/180.0);
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in deg2rad_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = d2r * X[n]; }

    return 0;
}


int deg2rad_d (double *Y, const double *X, const int N)
{
    const double d2r = M_PI / 180.0;
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in deg2rad_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = d2r * X[n]; }
    
    return 0;
}



int deg2rad_inplace_s (float *X, const int N)
{
    const float d2r = (float)(M_PI/180.0);

    //Checks
    if (N<0) { fprintf(stderr,"error in deg2rad_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    cblas_sscal(N,d2r,X,1);

    return 0;
}


int deg2rad_inplace_d (double *X, const int N)
{
    const double d2r = M_PI / 180.0;

    //Checks
    if (N<0) { fprintf(stderr,"error in deg2rad_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    cblas_dscal(N,d2r,X,1);
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
