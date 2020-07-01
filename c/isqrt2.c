//Sets all N elements of Y equal to isqrt2.
//For complex cases, only real part is set to isqrt2.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int isqrt2_s (float *Y, const int N);
int isqrt2_d (double *Y, const int N);
int isqrt2_c (float *Y, const int N);
int isqrt2_z (double *Y, const int N);


int isqrt2_s (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in isqrt2_s: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v = (float)M_SQRT1_2;

    cblas_scopy(N,&v,0,Y,1);

    return 0;
}


int isqrt2_d (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in isqrt2_d: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v = M_SQRT1_2;

    cblas_dcopy(N,&v,0,Y,1);

    return 0;
}


int isqrt2_c (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in isqrt2_c: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v[2] = {(float)M_SQRT1_2,0.0f};

    cblas_ccopy(N,v,0,Y,1);

    return 0;
}


int isqrt2_z (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in isqrt2_z: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v[2] = {M_SQRT1_2,0.0};

    cblas_zcopy(N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
