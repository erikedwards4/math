//Sets all N elements of Y equal to ipi.
//For complex cases, only real part is set to ipi.

#include <stdio.h>
//#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifndef M_1_PI
    #define M_1_PI 0.318309886183790671538
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ipi_s (float *Y, const size_t N);
int ipi_d (double *Y, const size_t N);
int ipi_c (float *Y, const size_t N);
int ipi_z (double *Y, const size_t N);


int ipi_s (float *Y, const size_t N)
{
    const float v = (float)M_1_PI;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int ipi_d (double *Y, const size_t N)
{
    const double v = M_1_PI;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int ipi_c (float *Y, const size_t N)
{
    const float v[2] = {(float)M_1_PI,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int ipi_z (double *Y, const size_t N)
{
    const double v[2] = {M_1_PI,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
