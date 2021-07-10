//Sets all N elements of Y equal to val.
//For complex cases, only real part is set to val.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int fill_s (float *Y, const size_t N, const float val);
int fill_d (double *Y, const size_t N, const double val);
int fill_c (float *Y, const size_t N, const float val);
int fill_z (double *Y, const size_t N, const double val);


int fill_s (float *Y, const size_t N, const float val)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = val; }

    return 0;
}


int fill_d (double *Y, const size_t N, const double val)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = val; }

    return 0;
}


int fill_c (float *Y, const size_t N, const float val)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = val; *++Y = 0.0f; }

    return 0;
}


int fill_z (double *Y, const size_t N, const double val)
{
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = val; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
