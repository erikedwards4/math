//Sets all N elements of Y equal to eps (FLT_EPSILON or DBL_EPSILON).
//For complex cases, only real part is set to eps.

#include <stdio.h>
#include <float.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int eps_s (float *Y, const size_t N);
int eps_d (double *Y, const size_t N);
int eps_c (float *Y, const size_t N);
int eps_z (double *Y, const size_t N);


int eps_s (float *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = FLT_EPSILON; }

    return 0;
}


int eps_d (double *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = DBL_EPSILON; }

    return 0;
}


int eps_c (float *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = FLT_EPSILON; *++Y = 0.0f; }

    return 0;
}


int eps_z (double *Y, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = DBL_EPSILON; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
