//Sets all N elements of Y equal to val.
//For complex cases, only real part is set to val.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int fill_s (float *Y, const size_t N, const float val);
int fill_d (double *Y, const size_t N, const double val);
int fill_c (float *Y, const size_t N, const float val);
int fill_z (double *Y, const size_t N, const double val);


int fill_s (float *Y, const size_t N, const float val)
{
    cblas_scopy((int)N,&val,0,Y,1);

    return 0;
}


int fill_d (double *Y, const size_t N, const double val)
{
    cblas_dcopy((int)N,&val,0,Y,1);

    return 0;
}


int fill_c (float *Y, const size_t N, const float val)
{
    const float v[2] = {val,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int fill_z (double *Y, const size_t N, const double val)
{
    const double v[2] = {val,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
