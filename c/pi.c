//Sets all N elements of Y equal to pi.
//For complex cases, only real part is set to pi.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int pi_s (float *Y, const int N);
int pi_d (double *Y, const int N);
int pi_c (float *Y, const int N);
int pi_z (double *Y, const int N);


int pi_s (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in pi_s: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v = (float)M_PI;

    cblas_scopy(N,&v,0,Y,1);

    return 0;
}


int pi_d (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in pi_d: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v = M_PI;

    cblas_dcopy(N,&v,0,Y,1);

    return 0;
}


int pi_c (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in pi_c: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v[2] = {(float)M_PI,0.0f};

    cblas_ccopy(N,v,0,Y,1);

    return 0;
}


int pi_z (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in pi_z: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v[2] = {M_PI,0.0};

    cblas_zcopy(N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
