//Sets all N elements of Y equal to 1/sqrt(2).
//For complex cases, only real part is set to isqrt2.

#include <stdio.h>
#include "codee_math.h"
//#include <math.h>

#ifndef M_SQRT1_2
    #define M_SQRT1_2 0.707106781186547524401
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int isqrt2_s (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_SQRT1_2; }

    return 0;
}


int isqrt2_d (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_SQRT1_2; }

    return 0;
}


int isqrt2_c (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_SQRT1_2; *++Y = 0.0f; }

    return 0;
}


int isqrt2_z (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_SQRT1_2; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
