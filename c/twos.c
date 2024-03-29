//Sets all N elements of Y equal to 2.
//For complex cases, only real part is set to 2.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int twos_s (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = 2.0f; }

    return 0;
}


int twos_d (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = 2.0; }

    return 0;
}


int twos_c (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = 2.0f; *++Y = 0.0f; }

    return 0;
}


int twos_z (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = 2.0; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
