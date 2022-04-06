//Sets all N elements of Y equal to ipi.
//For complex cases, only real part is set to ipi.

#include <stdio.h>
#include "codee_math.h"
//#include <math.h>

#ifndef M_1_PI
    #define M_1_PI 0.318309886183790671538
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int ipi_s (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_1_PI; }

    return 0;
}


int ipi_d (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_1_PI; }

    return 0;
}


int ipi_c (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_1_PI; *++Y = 0.0f; }

    return 0;
}


int ipi_z (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_1_PI; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
