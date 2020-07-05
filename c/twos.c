//Sets all N elements of Y equal to 2.
//For complex cases, only real part is set to 2.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int twos_s (float *Y, const size_t N);
int twos_d (double *Y, const size_t N);
int twos_c (float *Y, const size_t N);
int twos_z (double *Y, const size_t N);


int twos_s (float *Y, const size_t N)
{
    const float v = 2.0f;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    cblas_scopy((int)N,&v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int twos_d (double *Y, const size_t N)
{
    const double v = 2.0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    
    cblas_dcopy((int)N,&v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int twos_c (float *Y, const size_t N)
{
    const float v[2] = {2.0f,0.0f};
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    cblas_ccopy((int)N,v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int twos_z (double *Y, const size_t N)
{
    const double v[2] = {2.0,0.0};
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    cblas_zcopy((int)N,v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
