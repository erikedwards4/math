//Sets all N elements of Y equal to pi.
//For complex cases, only real part is set to pi.

#include <stdio.h>
#include "codee_math.h"
//#include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int pi_s (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_PI; }

    return 0;
}


int pi_d (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_PI; }

    return 0;
}


int pi_c (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_PI; *++Y = 0.0f; }

    return 0;
}


int pi_z (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_PI; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
