//Sets all N elements of Y equal to 1.
//For complex cases, only real part is set to 1.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ones_s (float *Y, const int N);
int ones_d (double *Y, const int N);
int ones_c (float *Y, const int N);
int ones_z (double *Y, const int N);


int ones_s (float *Y, const int N)
{
    const float o = 1.0f;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in ones_s: N (num elements Y) must be nonnegative\n"); return 1; }

    cblas_scopy(N,&o,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int ones_d (double *Y, const int N)
{
    const double o = 1.0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in ones_d: N (num elements Y) must be nonnegative\n"); return 1; }

    cblas_dcopy(N,&o,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int ones_c (float *Y, const int N)
{
    const float o[2] = {1.0f,0.0f};
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in ones_c: N (num elements Y) must be nonnegative\n"); return 1; }

    cblas_ccopy(N,o,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int ones_z (double *Y, const int N)
{
    const double o[2] = {1.0,0.0};
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in ones_z: N (num elements Y) must be nonnegative\n"); return 1; }

    cblas_zcopy(N,o,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
