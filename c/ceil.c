//Gets ceil of each element of X (rounds toward Inf).
//This has in-place and not-in-place versions.
//For complex input, ceils real and imag parts separately.

#include <stdio.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ceil_s (float *Y, const float *X, const size_t N);
int ceil_d (double *Y, const double *X, const size_t N);
int ceil_c (float *Y, const float *X, const size_t N);
int ceil_z (double *Y, const double *X, const size_t N);

int ceil_inplace_s (float *X, const size_t N);
int ceil_inplace_d (double *X, const size_t N);
int ceil_inplace_c (float *X, const size_t N);
int ceil_inplace_z (double *X, const size_t N);


int ceil_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = ceilf(*X); }

    return 0;
}


int ceil_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = ceil(*X); }
    
    return 0;
}


int ceil_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = ceilf(*X); }
    
    return 0;
}


int ceil_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = ceil(*X); }
    
    return 0;
}


int ceil_inplace_s (float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = ceilf(*X); }

    return 0;
}


int ceil_inplace_d (double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = ceil(*X); }
    
    return 0;
}


int ceil_inplace_c (float *X, const size_t N)
{
    for (size_t n=2u*N; n>0u; --n, ++X) { *X = ceilf(*X); }
    
    return 0;
}


int ceil_inplace_z (double *X, const size_t N)
{
    for (size_t n=2u*N; n>0u; --n, ++X) { *X = ceil(*X); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
