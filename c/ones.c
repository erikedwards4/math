//Sets all N elements of Y equal to 1.
//For complex cases, only real part is set to 1.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ones_s (float *Y, const size_t N);
int ones_d (double *Y, const size_t N);
int ones_c (float *Y, const size_t N);
int ones_z (double *Y, const size_t N);


int ones_s (float *Y, const size_t N)
{
    const float o = 1.0f;

    cblas_scopy((int)N,&o,0,Y,1);

    return 0;
}


int ones_d (double *Y, const size_t N)
{
    const double o = 1.0;

    cblas_dcopy((int)N,&o,0,Y,1);

    return 0;
}


int ones_c (float *Y, const size_t N)
{
    const float o[2] = {1.0f,0.0f};

    cblas_ccopy((int)N,o,0,Y,1);

    return 0;
}


int ones_z (double *Y, const size_t N)
{
    const double o[2] = {1.0,0.0};

    cblas_zcopy((int)N,o,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
