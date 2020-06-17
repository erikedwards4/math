//Sets all N elements of Y equal to val.
//For complex cases, only real part is set to val.

#include <stdio.h>
#include <float.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int fill_s (float *Y, const int N, const float val);
int fill_d (double *Y, const int N, const double val);
int fill_c (float *Y, const int N, const float val);
int fill_z (double *Y, const int N, const double val);


int fill_s (float *Y, const int N, const float val)
{
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in fill_s: N (num elements Y) must be nonnegative\n"); return 1; }

    cblas_scopy(N,&val,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int fill_d (double *Y, const int N, const double val)
{
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in fill_d: N (num elements Y) must be nonnegative\n"); return 1; }

    cblas_dcopy(N,&val,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int fill_c (float *Y, const int N, const float val)
{
    const float v[2] = {val,0.0f};
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in fill_c: N (num elements Y) must be nonnegative\n"); return 1; }

    cblas_ccopy(N,v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int fill_z (double *Y, const int N, const double val)
{
    const double v[2] = {val,0.0};
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in fill_z: N (num elements Y) must be nonnegative\n"); return 1; }

    cblas_zcopy(N,v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
