//Sets all N elements of Y equal to e.
//For complex cases, only real part is set to e.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int e_s (float *Y, const int N);
int e_d (double *Y, const int N);
int e_c (float *Y, const int N);
int e_z (double *Y, const int N);


int e_s (float *Y, const int N)
{
    const float v = (float)M_E;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in e_s: N (num elements Y) must be nonnegative\n"); return 1; }

    cblas_scopy(N,&v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int e_d (double *Y, const int N)
{
    const double v = M_E;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in e_d: N (num elements Y) must be nonnegative\n"); return 1; }

    cblas_dcopy(N,&v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int e_c (float *Y, const int N)
{
    const float v[2] = {(float)M_E,0.0f};
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in e_c: N (num elements Y) must be nonnegative\n"); return 1; }

    cblas_ccopy(N,v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int e_z (double *Y, const int N)
{
    const double v[2] = {M_E,0.0};
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<0) { fprintf(stderr,"error in e_z: N (num elements Y) must be nonnegative\n"); return 1; }

    cblas_zcopy(N,v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
