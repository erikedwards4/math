//Sets all N elements of Y equal to ln2 (M_LN2).
//For complex cases, only real part is set to ln2.

#include <stdio.h>
//#include <math.h>

#ifndef M_LN2
    #define M_LN2 0.693147180559945309417
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ln2_s (float *Y, const size_t N);
int ln2_d (double *Y, const size_t N);
int ln2_c (float *Y, const size_t N);
int ln2_z (double *Y, const size_t N);


int ln2_s (float *Y, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++Y) { *Y = (float)M_LN2; }

    return 0;
}


int ln2_d (double *Y, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++Y) { *Y = M_LN2; }

    return 0;
}


int ln2_c (float *Y, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++Y) { *Y = (float)M_LN2; *++Y = 0.0f; }

    return 0;
}


int ln2_z (double *Y, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++Y) { *Y = M_LN2; *++Y = 0.0; }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
