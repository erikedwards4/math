//Sets all N elements of Y equal to sqrt2.
//For complex cases, only real part is set to sqrt2.

#include <stdio.h>
//#include <math.h>

#ifndef M_SQRT2
    #define M_SQRT2 1.41421356237309504880
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sqrt2_s (float *Y, const size_t N);
int sqrt2_d (double *Y, const size_t N);
int sqrt2_c (float *Y, const size_t N);
int sqrt2_z (double *Y, const size_t N);


int sqrt2_s (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_SQRT2; }

    return 0;
}


int sqrt2_d (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_SQRT2; }

    return 0;
}


int sqrt2_c (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_SQRT2; *++Y = 0.0f; }

    return 0;
}


int sqrt2_z (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_SQRT2; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
