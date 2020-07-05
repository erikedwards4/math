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

int nan_s (float *Y, const size_t N);
int nan_d (double *Y, const size_t N);
int nan_c (float *Y, const size_t N);
int nan_z (double *Y, const size_t N);


int nan_s (float *Y, const size_t N)
{
    const float v = NAN;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int nan_d (double *Y, const size_t N)
{
    const double v = (double)NAN;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int nan_c (float *Y, const size_t N)
{
    const float v[2] = {NAN,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int nan_z (double *Y, const size_t N)
{
    const double v[2] = {(double)NAN,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
