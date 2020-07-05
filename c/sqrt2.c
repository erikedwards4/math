//Sets all N elements of Y equal to sqrt2.
//For complex cases, only real part is set to sqrt2.

#include <stdio.h>
//#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifndef M_SQRT2
    #define M_SQRT2 1.41421356237309504880
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sqrt2_s (float *Y, const size_t N);
int sqrt2_d (double *Y, const size_t N);
int sqrt2_c (float *Y, const size_t N);
int sqrt2_z (double *Y, const size_t N);


int sqrt2_s (float *Y, const size_t N)
{
    const float v = (float)M_SQRT2;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    cblas_scopy((int)N,&v,0,Y,1);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int sqrt2_d (double *Y, const size_t N)
{
    const double v = M_SQRT2;

    cblas_dcopy((int)N,&v,0,Y,1);

    return 0;
}


int sqrt2_c (float *Y, const size_t N)
{
    const float v[2] = {(float)M_SQRT2,0.0f};

    cblas_ccopy((int)N,v,0,Y,1);

    return 0;
}


int sqrt2_z (double *Y, const size_t N)
{
    const double v[2] = {M_SQRT2,0.0};

    cblas_zcopy((int)N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
