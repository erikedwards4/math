//Sets all N elements of Y equal to pi/2.
//For complex cases, only real part is set to pi/2.

#include <stdio.h>
#include "codee_math.h"
//#include <math.h>

#ifndef M_PI_2
    #define M_PI_2 1.57079632679489661923
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int pi_2_s (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_PI_2; }

    return 0;
}


int pi_2_d (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_PI_2; }

    return 0;
}


int pi_2_c (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_PI_2; *++Y = 0.0f; }

    return 0;
}


int pi_2_z (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_PI_2; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
