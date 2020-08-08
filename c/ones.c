//Sets all N elements of Y equal to 1.
//For complex cases, only real part is set to 1.

//The for loop is much faster for small N,
//and slightly faster for large N, compared to the cblas_?copy.
//E.g., for N=1e8, cblas_scopy takes ~103 ms,
//and for loope takes ~103 ms at O1 level of optimization;
//but for loop takes ~83 ms at O2 and O3 levels of optimization.

#include <stdio.h>
//#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ones_s (float *Y, const size_t N);
int ones_d (double *Y, const size_t N);
int ones_c (float *Y, const size_t N);
int ones_z (double *Y, const size_t N);


int ones_s (float *Y, const size_t N)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    //const float o = 1.0f;
    //cblas_scopy((int)N,&o,0,Y,1);

    for (size_t n=0; n<N; ++n, ++Y) { *Y = 1.0f; }

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int ones_d (double *Y, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++Y) { *Y = 1.0; }

    return 0;
}


int ones_c (float *Y, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++Y) { *Y = 1.0f; *++Y = 0.0f; }

    return 0;
}


int ones_z (double *Y, const size_t N)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    for (size_t n=0; n<N; ++n, ++Y) { *Y = 1.0; *++Y = 0.0; }

    //const double o[2] = {1.0,0.0};
    //cblas_zcopy((int)N,o,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
