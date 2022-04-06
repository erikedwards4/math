//Sets all N elements of Y equal to realmin (FLT_MIN or DBL_MIN).
//This is a small, positive number (smaller than FLT_EPSILON or DBL_EPSILON).
//For complex cases, only real part is set to realmin.

#include <stdio.h>
#include <float.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int realmin_s (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = FLT_MIN; }

    return 0;
}


int realmin_d (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = DBL_MIN; }

    return 0;
}


int realmin_c (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = FLT_MIN; *++Y = 0.0f; }

    return 0;
}


int realmin_z (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = DBL_MIN; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
