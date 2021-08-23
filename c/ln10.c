//Sets all N elements of Y equal to ln10 (M_LN10).
//For complex cases, only real part is set to ln10.

#include <stdio.h>
//#include <math.h>

#ifndef M_LN10
    #define M_LN10 2.30258509299404568402
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ln10_s (float *Y, const size_t N);
int ln10_d (double *Y, const size_t N);
int ln10_c (float *Y, const size_t N);
int ln10_z (double *Y, const size_t N);


int ln10_s (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_LN10; }

    return 0;
}


int ln10_d (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_LN10; }

    return 0;
}


int ln10_c (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_LN10; *++Y = 0.0f; }

    return 0;
}


int ln10_z (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_LN10; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
