//Gets real part of complex-valued input X.
//This has in-place and not-in-place versions.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int real_c (float *Y, const float *X, const size_t N);
int real_z (double *Y, const double *X, const size_t N);


int real_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, X+=2, ++Y) { *Y = *X; }

    return 0;
}


int real_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, X+=2, ++Y) { *Y = *X; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
