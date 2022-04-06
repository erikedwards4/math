//Gets floor of each element of X (rounds toward -Inf).
//This has in-place and not-in-place versions.
//For complex input, floors real and imag parts separately.

#include <stdio.h>
#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int floor_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = floorf(*X); }

    return 0;
}


int floor_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = floor(*X); }
    
    return 0;
}


int floor_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = floorf(*X); }
    
    return 0;
}


int floor_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = floor(*X); }
    
    return 0;
}


int floor_inplace_s (float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = floorf(*X); }

    return 0;
}


int floor_inplace_d (double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = floor(*X); }
    
    return 0;
}


int floor_inplace_c (float *X, const size_t N)
{
    for (size_t n=2u*N; n>0u; --n, ++X) { *X = floorf(*X); }
    
    return 0;
}


int floor_inplace_z (double *X, const size_t N)
{
    for (size_t n=2u*N; n>0u; --n, ++X) { *X = floor(*X); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
