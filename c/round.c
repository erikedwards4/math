//Rounds of each element of X (towards nearest integer).
//This has in-place and not-in-place versions.
//For complex input, rounds real and imag parts separately.

#include <stdio.h>
#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int round_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = roundf(*X); }

    return 0;
}


int round_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = round(*X); }
    
    return 0;
}


int round_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = roundf(*X); }
    
    return 0;
}


int round_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = round(*X); }
    
    return 0;
}


int round_inplace_s (float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = roundf(*X); }

    return 0;
}


int round_inplace_d (double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = round(*X); }
    
    return 0;
}


int round_inplace_c (float *X, const size_t N)
{
    for (size_t n=2u*N; n>0u; --n, ++X) { *X = roundf(*X); }
    
    return 0;
}


int round_inplace_z (double *X, const size_t N)
{
    for (size_t n=2u*N; n>0u; --n, ++X) { *X = round(*X); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
