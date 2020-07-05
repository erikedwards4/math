//Sets all N elements of Y equal to isqrt2.
//For complex cases, only real part is set to isqrt2.

#include <stdio.h>
//#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifndef M_SQRT1_2
    #define M_SQRT1_2 0.707106781186547524401
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int isqrt2_s (float *Y, const size_t N);
int isqrt2_d (double *Y, const size_t N);
int isqrt2_c (float *Y, const size_t N);
int isqrt2_z (double *Y, const size_t N);


int isqrt2_s (float *Y, const size_t N)
{
    const float v = (float)M_SQRT1_2;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int isqrt2_d (double *Y, const size_t N)
{
    const double v = M_SQRT1_2;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int isqrt2_c (float *Y, const size_t N)
{
    const float v[2] = {(float)M_SQRT1_2,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int isqrt2_z (double *Y, const size_t N)
{
    const double v[2] = {M_SQRT1_2,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
