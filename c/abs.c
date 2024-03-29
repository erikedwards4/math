//This just gets absolute value of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int abs_s (float *Y, const float *X, const size_t N)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&tic);
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = fabsf(*X); }
    //for (size_t n=N; n>0u; --n) { Y[n] = fabsf(X[n]); }                 //same speed, but more instructions
    //for (size_t n=N; n>0u; --n) { Y[n] = (X[n]<0.0f) ? -X[n] : X[n]; }  //slower
    //clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int abs_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = fabs(*X); }
    
    return 0;
}


int abs_c (float *Y, const float *X, const size_t N)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    for (size_t n=N; n>0u; --n, X+=2, ++Y)
    {
        *Y = sqrtf(*X**X + *(X+1)**(X+1));
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int abs_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, X+=2, ++Y)
    {
        *Y = sqrt(*X**X + *(X+1)**(X+1));
    }
    
    return 0;
}


int abs_inplace_s (float *X, const size_t N)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    for (size_t n=N; n>0u; --n, ++X) { *X = fabsf(*X); }
    //for (size_t n=N; n>0u; --n) { if (X[n]<0.0f) { X[n] = -X[n]; } }
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int abs_inplace_d (double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = fabs(*X); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
