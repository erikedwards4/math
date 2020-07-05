//Sets all N elements of Y equal to e.
//For complex cases, only real part is set to e.

#include <stdio.h>
//#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifndef M_E
    #define M_E 2.71828182845904523536
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int e_s (float *Y, const size_t N);
int e_d (double *Y, const size_t N);
int e_c (float *Y, const size_t N);
int e_z (double *Y, const size_t N);


int e_s (float *Y, const size_t N)
{
    const float v = (float)M_E;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int e_d (double *Y, const size_t N)
{
    const double v = M_E;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int e_c (float *Y, const size_t N)
{
    const float v[2] = {(float)M_E,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int e_z (double *Y, const size_t N)
{
    const double v[2] = {M_E,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
