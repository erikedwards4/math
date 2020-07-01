//Sets all N elements of Y equal to nan.
//For complex cases, only real part is set to nan (NAN).

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int nan_s (float *Y, const int N);
int nan_d (double *Y, const int N);
int nan_c (float *Y, const int N);
int nan_z (double *Y, const int N);


int nan_s (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in nan_s: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v = NAN;

    cblas_scopy(N,&v,0,Y,1);

    return 0;
}


int nan_d (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in nan_d: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v = (double)NAN;

    cblas_dcopy(N,&v,0,Y,1);

    return 0;
}


int nan_c (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in nan_c: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v[2] = {NAN,0.0f};

    cblas_ccopy(N,v,0,Y,1);

    return 0;
}


int nan_z (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in nan_z: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v[2] = {(double)NAN,0.0};

    cblas_zcopy(N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
