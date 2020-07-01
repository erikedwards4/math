//Sets all N elements of Y equal to 1.
//For complex cases, only real part is set to 1.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ones_s (float *Y, const int N);
int ones_d (double *Y, const int N);
int ones_c (float *Y, const int N);
int ones_z (double *Y, const int N);


int ones_s (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in ones_s: N (num elements Y) must be nonnegative\n"); return 1; }

    const float o = 1.0f;

    cblas_scopy(N,&o,0,Y,1);

    return 0;
}


int ones_d (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in ones_d: N (num elements Y) must be nonnegative\n"); return 1; }

    const double o = 1.0;

    cblas_dcopy(N,&o,0,Y,1);

    return 0;
}


int ones_c (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in ones_c: N (num elements Y) must be nonnegative\n"); return 1; }

    const float o[2] = {1.0f,0.0f};

    cblas_ccopy(N,o,0,Y,1);

    return 0;
}


int ones_z (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in ones_z: N (num elements Y) must be nonnegative\n"); return 1; }

    const double o[2] = {1.0,0.0};

    cblas_zcopy(N,o,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
