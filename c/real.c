//Gets real part of complex-valued input X.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int real_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }

    return 0;
}


int real_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }

    return 0;
}


int real_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, X+=2, ++Y) { *Y = *X; }

    return 0;
}


int real_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, X+=2, ++Y) { *Y = *X; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
