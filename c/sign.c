//This gets sign for each element of X.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int sign_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = (float)((*X>0.0f)-(*X<0.0f)); }

    return 0;
}


int sign_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = (double)((*X>0.0)-(*X<0.0)); }
    
    return 0;
}


int sign_inplace_s (float *X, const size_t N)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    for (size_t n=N; n>0u; --n, ++X) { *X = (float)((*X>0.0f)-(*X<0.0f)); }
    // for (size_t n=N; n>0u; --n, ++X)
    // {
    //     if (*X<0.0f) { *X = -1.0f; }
    //     else if (*X>0.0f) { *X = 1.0f; }
    //     else { *X = 0.0f; }
    // }

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int sign_inplace_d (double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = (double)((*X>0.0)-(*X<0.0)); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
