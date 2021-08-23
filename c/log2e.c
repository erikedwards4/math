//Sets all N elements of Y equal to log2e (M_LOG2E).
//For complex cases, only real part is set to log2e.

#include <stdio.h>
//#include <math.h>

#ifndef M_LOG2E
    #define M_LOG2E 1.44269504088896340736
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int log2e_s (float *Y, const size_t N);
int log2e_d (double *Y, const size_t N);
int log2e_c (float *Y, const size_t N);
int log2e_z (double *Y, const size_t N);


int log2e_s (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_LOG2E; }

    return 0;
}


int log2e_d (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_LOG2E; }

    return 0;
}


int log2e_c (float *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = (float)M_LOG2E; *++Y = 0.0f; }

    return 0;
}


int log2e_z (double *Y, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++Y) { *Y = M_LOG2E; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
