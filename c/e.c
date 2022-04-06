//Sets all N elements of Y equal to e.
//For complex cases, only real part is set to e.

#include <stdio.h>
#include "codee_math.h"
//#include <math.h>

#ifndef M_E
    #define M_E 2.71828182845904523536
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int e_s (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_E; }

    return 0;
}


int e_d (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_E; }

    return 0;
}


int e_c (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_E; *++Y = 0.0f; }

    return 0;
}


int e_z (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_E; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
