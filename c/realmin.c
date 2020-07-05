//Sets all N elements of Y equal to realmin (FLT_MIN or DBL_MIN).
//This is a small, positive number (smaller than FLT_EPSILON or DBL_EPSILON).
//For complex cases, only real part is set to realmin.

#include <stdio.h>
#include <float.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int realmin_s (float *Y, const size_t N);
int realmin_d (double *Y, const size_t N);
int realmin_c (float *Y, const size_t N);
int realmin_z (double *Y, const size_t N);


int realmin_s (float *Y, const size_t N)
{
    const float v = FLT_MIN;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int realmin_d (double *Y, const size_t N)
{
    const double v = DBL_MIN;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int realmin_c (float *Y, const size_t N)
{
    const float v[2] = {FLT_MIN,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int realmin_z (double *Y, const size_t N)
{
    const double v[2] = {DBL_MIN,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
