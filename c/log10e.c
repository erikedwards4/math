//Sets all N elements of Y equal to log10e (M_LOG10E).
//For complex cases, only real part is set to log10e.

#include <stdio.h>
//#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifndef M_LOG10E
    #define M_LOG10E 0.434294481903251827651
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int log10e_s (float *Y, const size_t N);
int log10e_d (double *Y, const size_t N);
int log10e_c (float *Y, const size_t N);
int log10e_z (double *Y, const size_t N);


int log10e_s (float *Y, const size_t N)
{
    const float v = (float)M_LOG10E;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int log10e_d (double *Y, const size_t N)
{
    const double v = M_LOG10E;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int log10e_c (float *Y, const size_t N)
{
    const float v[2] = {(float)M_LOG10E,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int log10e_z (double *Y, const size_t N)
{
    const double v[2] = {M_LOG10E,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
