//Sets all N elements of Y equal to ln10 (M_LN10).
//For complex cases, only real part is set to ln10.

#include <stdio.h>
//#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifndef M_LN10
    #define M_LN10 2.30258509299404568402
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ln10_s (float *Y, const size_t N);
int ln10_d (double *Y, const size_t N);
int ln10_c (float *Y, const size_t N);
int ln10_z (double *Y, const size_t N);


int ln10_s (float *Y, const size_t N)
{
    const float v = (float)M_LN10;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int ln10_d (double *Y, const size_t N)
{
    const double v = M_LN10;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int ln10_c (float *Y, const size_t N)
{
    const float v[2] = {(float)M_LN10,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int ln10_z (double *Y, const size_t N)
{
    const double v[2] = {M_LN10,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
