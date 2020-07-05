//Sets all N elements of Y equal to ln2 (M_LN2).
//For complex cases, only real part is set to ln2.

#include <stdio.h>
//#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifndef M_LN2
    #define M_LN2 0.693147180559945309417
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ln2_s (float *Y, const size_t N);
int ln2_d (double *Y, const size_t N);
int ln2_c (float *Y, const size_t N);
int ln2_z (double *Y, const size_t N);


int ln2_s (float *Y, const size_t N)
{
    const float v = (float)M_LN2;

    cblas_scopy((int)N,&v,0,Y,1);

    return 0;
}


int ln2_d (double *Y, const size_t N)
{
    const double v = M_LN2;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int ln2_c (float *Y, const size_t N)
{
    const float v[2] = {(float)M_LN2,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int ln2_z (double *Y, const size_t N)
{
    const double v[2] = {M_LN2,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
