//For each element of X, does x++ (increment by 1).
//Since the actual operator is meant for integers, this does x += 1.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int plusplus_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X + 1.0f; }

    return 0;
}


int plusplus_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X + 1.0; }
    
    return 0;
}


int plusplus_inplace_s (float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X += 1.0f; }

    return 0;
}


int plusplus_inplace_d (double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X += 1.0; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
