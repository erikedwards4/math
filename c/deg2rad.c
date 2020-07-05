//Converts degrees to radians for each element of X.
//This has in-place and not-in-place versions.

//LAPACKE_slascl was definitely slower than cblas_sscal.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <lapacke.h>
//#include <time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int deg2rad_s (float *Y, const float *X, const size_t N);
int deg2rad_d (double *Y, const double *X, const size_t N);

int deg2rad_inplace_s (float *X, const size_t N);
int deg2rad_inplace_d (double *X, const size_t N);


int deg2rad_s (float *Y, const float *X, const size_t N)
{
    const float d2r = (float)(M_PI/180.0);

    for (size_t n=0; n<N; n++) { Y[n] = d2r * X[n]; }

    return 0;
}


int deg2rad_d (double *Y, const double *X, const size_t N)
{
    const double d2r = M_PI / 180.0;

    for (size_t n=0; n<N; n++) { Y[n] = d2r * X[n]; }
    
    return 0;
}



int deg2rad_inplace_s (float *X, const size_t N)
{
    const float d2r = (float)(M_PI/180.0);

    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    cblas_sscal((int)N,d2r,X,1);
    //if (LAPACKE_slascl_work(LAPACK_COL_MAJOR,'G',0,0,180.0f,(float)M_PI,1,N,X,1))
    //{ fprintf(stderr,"error in deg2rad_inplace_s: problem with LAPACKE function\n"); return 1; }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int deg2rad_inplace_d (double *X, const size_t N)
{
    const double d2r = M_PI / 180.0;

    cblas_dscal((int)N,d2r,X,1);
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
