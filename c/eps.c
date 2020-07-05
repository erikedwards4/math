//Sets all N elements of Y equal to eps (FLT_EPSILON or DBL_EPSILON).
//For complex cases, only real part is set to eps.

#include <stdio.h>
#include <float.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int eps_s (float *Y, const size_t N);
int eps_d (double *Y, const size_t N);
int eps_c (float *Y, const size_t N);
int eps_z (double *Y, const size_t N);


int eps_s (float *Y, const size_t N)
{
    const float v = FLT_EPSILON;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int eps_d (double *Y, const size_t N)
{
    const double v = DBL_EPSILON;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int eps_c (float *Y, const size_t N)
{
    const float v[2] = {FLT_EPSILON,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int eps_z (double *Y, const size_t N)
{
    const double v[2] = {DBL_EPSILON,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
