//Sets all N elements of Y equal to nan.
//For complex cases, only real part is set to nan (NAN).

#include <stdio.h>
#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int nan_s (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = NAN; }

    return 0;
}


int nan_d (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (double)NAN; }

    return 0;
}


int nan_c (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = NAN; *++Y = 0.0f; }

    return 0;
}


int nan_z (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (double)NAN; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
