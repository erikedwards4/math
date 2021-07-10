//Sets all N elements of Y equal to 0.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int zeros_s (float *Y, const size_t N);
int zeros_d (double *Y, const size_t N);
int zeros_c (float *Y, const size_t N);
int zeros_z (double *Y, const size_t N);


int zeros_s (float *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }

    return 0;
}


int zeros_d (double *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }

    return 0;
}


int zeros_c (float *Y, const size_t N)
{
    for (size_t n=0u; n<2*N; ++n, ++Y) { *Y = 0.0f; }

    return 0;
}


int zeros_z (double *Y, const size_t N)
{
    for (size_t n=0u; n<2*N; ++n, ++Y) { *Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
