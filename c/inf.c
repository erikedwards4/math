//Sets all N elements of Y equal to infinity (HUGE_VAL).
//For complex cases, only real part is set to infinity (HUGE_VAL).

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int inf_s (float *Y, const size_t N);
int inf_d (double *Y, const size_t N);
int inf_c (float *Y, const size_t N);
int inf_z (double *Y, const size_t N);


int inf_s (float *Y, const size_t N)
{
    const float v = HUGE_VALF;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int inf_d (double *Y, const size_t N)
{
    const double v = HUGE_VAL;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int inf_c (float *Y, const size_t N)
{
    const float v[2] = {HUGE_VALF,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int inf_z (double *Y, const size_t N)
{
    const double v[2] = {HUGE_VAL,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
