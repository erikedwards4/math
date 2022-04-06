//Sets all N elements of Y equal to pi/4.
//For complex cases, only real part is set to pi/4.

#include <stdio.h>
#include "codee_math.h"
//#include <math.h>

#ifndef M_PI_4
    #define M_PI_4 0.785398163397448309616
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int pi_4_s (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_PI_4; }

    return 0;
}


int pi_4_d (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_PI_4; }

    return 0;
}


int pi_4_c (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_PI_4; *++Y = 0.0f; }

    return 0;
}


int pi_4_z (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_PI_4; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
