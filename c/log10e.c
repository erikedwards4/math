//Sets all N elements of Y equal to log10e (M_LOG10E).
//For complex cases, only real part is set to log10e.

#include <stdio.h>
//#include <math.h>

#ifndef M_LOG10E
    #define M_LOG10E 0.434294481903251827651
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int log10e_s (float *Y, const size_t N);
int log10e_d (double *Y, const size_t N);
int log10e_c (float *Y, const size_t N);
int log10e_z (double *Y, const size_t N);


int log10e_s (float *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = (float)M_LOG10E; }

    return 0;
}


int log10e_d (double *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = M_LOG10E; }

    return 0;
}


int log10e_c (float *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = (float)M_LOG10E; *++Y = 0.0f; }

    return 0;
}


int log10e_z (double *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = M_LOG10E; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
