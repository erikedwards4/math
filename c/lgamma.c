//Gets log-gamma-function of each element of X.
//This has in-place and not-in-place versions.
//For complex input, gets lgamma of real and imag parts separately.

#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace openn {
extern "C" {
#endif

int lgamma_s (float *Y, const float *X, const int N);
int lgamma_d (double *Y, const double *X, const int N);
int lgamma_c (float *Y, const float *X, const int N);
int lgamma_z (double *Y, const double *X, const int N);

int lgamma_inplace_s (float *X, const int N);
int lgamma_inplace_d (double *X, const int N);
int lgamma_inplace_c (float *X, const int N);
int lgamma_inplace_z (double *X, const int N);


int lgamma_s (float *Y, const float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in lgamma_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = lgammaf(X[n]); }

    return 0;
}


int lgamma_d (double *Y, const double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in lgamma_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = lgamma(X[n]); }
    
    return 0;
}


int lgamma_c (float *Y, const float *X, const int N)
{
    int n;
    //struct timespec tic, toc;

    //Checks
    if (N<0) { fprintf(stderr,"error in lgamma_c: N (num elements X) must be nonnegative\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);
    for (n=0; n<2*N; n++) { Y[n] = lgammaf(X[n]); }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int lgamma_z (double *Y, const double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in lgamma_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<2*N; n++) { Y[n] = lgamma(X[n]); }
    
    return 0;
}


int lgamma_inplace_s (float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in lgamma_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { X[n] = lgammaf(X[n]); }

    return 0;
}


int lgamma_inplace_d (double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in lgamma_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { X[n] = lgamma(X[n]); }
    
    return 0;
}


int lgamma_inplace_c (float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in lgamma_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<2*N; n++) { X[n] = lgammaf(X[n]); }
    
    return 0;
}


int lgamma_inplace_z (double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in lgamma_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<2*N; n++) { X[n] = lgamma(X[n]); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
