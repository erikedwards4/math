//Gets real part of complex-valued input X.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int real_c (float *Y, const float *X, const size_t N);
int real_z (double *Y, const double *X, const size_t N);

int real_inplace_c (float *X, const size_t N);
int real_inplace_z (double *X, const size_t N);


int real_c (float *Y, const float *X, const size_t N)
{
    cblas_scopy((int)N,X,2,Y,1);

    return 0;
}


int real_z (double *Y, const double *X, const size_t N)
{
    cblas_dcopy((int)N,X,2,Y,1);
    
    return 0;
}


int real_inplace_c (float *X, const size_t N)
{
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    //cblas_scopy(N-1,&X[2],2,&X[1],1);
    for (size_t n=1, n2=2; n<N; ++n, n2+=2) { X[n] = X[n2]; }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int real_inplace_z (double *X, const size_t N)
{
    for (size_t n=1, n2=2; n<N; ++n, n2+=2) { X[n] = X[n2]; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
