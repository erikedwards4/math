//Sets all N elements of Y equal to 0.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int zeros_s (float *Y, const int N);
int zeros_d (double *Y, const int N);
int zeros_c (float *Y, const int N);
int zeros_z (double *Y, const int N);


int zeros_s (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in zeros_s: N (num elements Y) must be nonnegative\n"); return 1; }

    const float z = 0.0f;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    
    cblas_scopy(N,&z,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int zeros_d (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in zeros_d: N (num elements Y) must be nonnegative\n"); return 1; }

    const double z = 0.0;
    
    cblas_dcopy(N,&z,0,Y,1);

    return 0;
}


int zeros_c (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in zeros_c: N (num elements Y) must be nonnegative\n"); return 1; }

    const float z[2] = {0.0f,0.0f};

    cblas_ccopy(N,z,0,Y,1);

    return 0;
}


int zeros_z (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in zeros_z: N (num elements Y) must be nonnegative\n"); return 1; }

    const double z[2] = {0.0,0.0};

    cblas_zcopy(N,z,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
