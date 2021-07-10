//Sets all N elements of Y equal to realmax (FLT_MAX or DBL_MAX).
//For complex cases, only real part is set to realmax.

#include <stdio.h>
#include <float.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int realmax_s (float *Y, const size_t N);
int realmax_d (double *Y, const size_t N);
int realmax_c (float *Y, const size_t N);
int realmax_z (double *Y, const size_t N);


int realmax_s (float *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = FLT_MAX; }

    return 0;
}


int realmax_d (double *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = DBL_MAX; }

    return 0;
}


int realmax_c (float *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = FLT_MAX; *++Y = 0.0f; }

    return 0;
}


int realmax_z (double *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = DBL_MAX; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
