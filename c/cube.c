//This just cubes input X element-wise.
//For complex data, this is |X|.^3.

//This has in-place and not-in-place versions.

//powf(*X,3.0f) tested at same speed (perhaps a tiny bit slower?) as *X**X**X
//powf(*X,3.0f) has more assembly instructions at -O1 and -O2, but fewer at -O3.

#include <stdio.h>
#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int cube_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * *X * *X; }

    return 0;
}


int cube_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * *X * *X; }
    
    return 0;
}


int cube_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, X+=2, ++Y) { *Y = *X**X + *(X+1)**(X+1); *Y *= sqrtf(*Y); }
    
    return 0;
}


int cube_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, X+=2, ++Y) { *Y = *X**X + *(X+1)**(X+1); *Y *= sqrt(*Y); }
    
    return 0;
}


int cube_inplace_s (float *X, const size_t N)
{
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    for (size_t n=N; n>0u; --n, ++X) { *X = *X * *X * *X; }
    //for (size_t n=N; n>0u; --n, ++X) { *X = powf(*X,3.0f); } //same speed, but has few assembly lines at -O3
    //for (size_t n=N; n>0u; --n, ++X) { *X *= *X * *X; }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int cube_inplace_d (double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = *X * *X * *X; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
