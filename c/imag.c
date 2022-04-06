//Gets imaginary part of complex-valued input X (output is real-valued).

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int imag_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *++X; }

    return 0;
}


int imag_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *++X; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
