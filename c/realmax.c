//Sets all N elements of Y equal to realmax (FLT_MAX or DBL_MAX).
//For complex cases, only real part is set to realmax.

#include <stdio.h>
#include <float.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int realmax_s (float *Y, const size_t N);
int realmax_d (double *Y, const size_t N);
int realmax_c (float *Y, const size_t N);
int realmax_z (double *Y, const size_t N);


int realmax_s (float *Y, const size_t N)
{
    const float v = FLT_MAX;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int realmax_d (double *Y, const size_t N)
{
    const double v = DBL_MAX;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int realmax_c (float *Y, const size_t N)
{
    const float v[2] = {FLT_MAX,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int realmax_z (double *Y, const size_t N)
{
    const double v[2] = {DBL_MAX,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
