//Sets all N elements of Y equal to infinity (HUGE_VAL).
//For complex cases, only real part is set to infinity (HUGE_VAL).

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int inf_s (float *Y, const size_t N);
int inf_d (double *Y, const size_t N);
int inf_c (float *Y, const size_t N);
int inf_z (double *Y, const size_t N);


int inf_s (float *Y, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++Y) { *Y = HUGE_VALF; }

    return 0;
}


int inf_d (double *Y, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++Y) { *Y = HUGE_VAL; }

    return 0;
}


int inf_c (float *Y, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++Y) { *Y = HUGE_VALF; *++Y = 0.0f; }

    return 0;
}


int inf_z (double *Y, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++Y) { *Y = HUGE_VAL; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
