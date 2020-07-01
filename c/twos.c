//Sets all N elements of Y equal to 2.
//For complex cases, only real part is set to 2.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int twos_s (float *Y, const int N);
int twos_d (double *Y, const int N);
int twos_c (float *Y, const int N);
int twos_z (double *Y, const int N);


int twos_s (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in twos_s: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v = 2.0f;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    cblas_scopy(N,&v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int twos_d (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in twos_d: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v = 2.0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    
    cblas_dcopy(N,&v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int twos_c (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in twos_c: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v[2] = {2.0f,0.0f};
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    cblas_ccopy(N,v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int twos_z (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in twos_z: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v[2] = {2.0,0.0};
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    cblas_zcopy(N,v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
