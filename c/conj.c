//Gets conj part of complex-valued input X.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int conj_c (float *Y, const float *X, const int N);
int conj_z (double *Y, const double *X, const int N);

int conj_inplace_c (float *X, const int N);
int conj_inplace_z (double *X, const int N);


int conj_c (float *Y, const float *X, const int N)
{
    //Checks
    if (N<0) { fprintf(stderr,"error in conj_s: N (num elements X) must be nonnegative\n"); return 1; }

    cblas_scopy(2*N,X,1,Y,1);
    cblas_sscal(N,-1.0f,&Y[1],2);

    return 0;
}


int conj_z (double *Y, const double *X, const int N)
{
    //struct timespec tic, toc;

    //Checks
    if (N<0) { fprintf(stderr,"error in conj_z: N (num elements X) must be nonnegative\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);
    cblas_dcopy(2*N,X,1,Y,1);
    cblas_dscal(N,-1.0,&Y[1],2);
    //for (n=0; n<2*N; n+=2) { Y[n] = X[n]; }
    //for (n=1; n<2*N; n+=2) { Y[n] = -X[n]; }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int conj_inplace_c (float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in conj_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=1; n<N; n+=2) { X[n] = -X[n]; }
    
    return 0;
}


int conj_inplace_z (double *X, const int N)
{
    int n;
    //struct timespec tic, toc;

    //Checks
    if (N<0) { fprintf(stderr,"error in conj_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);
    //cblas_dscal(N,-1.0,&X[1],2);
    for (n=1; n<N; n+=2) { X[n] = -X[n]; }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
