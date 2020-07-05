//Sets all N elements of Y equal to log2e (M_LOG2E).
//For complex cases, only real part is set to log2e.

#include <stdio.h>
//#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifndef M_LOG2E
    #define M_LOG2E 1.44269504088896340736
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int log2e_s (float *Y, const size_t N);
int log2e_d (double *Y, const size_t N);
int log2e_c (float *Y, const size_t N);
int log2e_z (double *Y, const size_t N);


int log2e_s (float *Y, const size_t N)
{
    const float v = (float)M_LOG2E;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int log2e_d (double *Y, const size_t N)
{
    const double v = M_LOG2E;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int log2e_c (float *Y, const size_t N)
{
    const float v[2] = {(float)M_LOG2E,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int log2e_z (double *Y, const size_t N)
{
    const double v[2] = {M_LOG2E,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
