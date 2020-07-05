//Sets all N elements of Y equal to pi.
//For complex cases, only real part is set to pi.

#include <stdio.h>
//#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int pi_s (float *Y, const size_t N);
int pi_d (double *Y, const size_t N);
int pi_c (float *Y, const size_t N);
int pi_z (double *Y, const size_t N);


int pi_s (float *Y, const size_t N)
{
    const float v = (float)M_PI;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int pi_d (double *Y, const size_t N)
{
    const double v = M_PI;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int pi_c (float *Y, const size_t N)
{
    const float v[2] = {(float)M_PI,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int pi_z (double *Y, const size_t N)
{
    const double v[2] = {M_PI,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
