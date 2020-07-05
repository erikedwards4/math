//Sets all N elements of Y equal to pi/2.
//For complex cases, only real part is set to pi/2.

#include <stdio.h>
//#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifndef M_PI_2
    #define M_PI_2 1.57079632679489661923
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int pi_2_s (float *Y, const size_t N);
int pi_2_d (double *Y, const size_t N);
int pi_2_c (float *Y, const size_t N);
int pi_2_z (double *Y, const size_t N);


int pi_2_s (float *Y, const size_t N)
{
    const float v = (float)M_PI_2;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int pi_2_d (double *Y, const size_t N)
{
    const double v = M_PI_2;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int pi_2_c (float *Y, const size_t N)
{
    const float v[2] = {(float)M_PI_2,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int pi_2_z (double *Y, const size_t N)
{
    const double v[2] = {M_PI_2,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
