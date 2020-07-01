//Gets proj part of complex-valued input X.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <complex.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int proj_c (float *Y, const float *X, const int N);
int proj_z (double *Y, const double *X, const int N);

int proj_inplace_c (float *X, const int N);
int proj_inplace_z (double *X, const int N);


int proj_c (float *Y, const float *X, const int N)
{
    int n;
    _Complex float y;

    //Checks
    if (N<0) { fprintf(stderr,"error in proj_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<2*N; n+=2)
    {
        y = cprojf(X[n]+1.0if*X[n+1]);
        memcpy(&Y[n],(float *)&y,2*sizeof(float));
    }

    return 0;
}


int proj_z (double *Y, const double *X, const int N)
{
    int n;
    _Complex double y;

    //Checks
    if (N<0) { fprintf(stderr,"error in proj_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<2*N; n+=2)
    {
        y = cproj(X[n]+1.0i*X[n+1]);
        memcpy(&Y[n],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int proj_inplace_c (float *X, const int N)
{
    int n;
    _Complex float y;

    //Checks
    if (N<0) { fprintf(stderr,"error in proj_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<2*N; n+=2)
    {
        y = cprojf(X[n]+1.0if*X[n+1]);
        memcpy(&X[n],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int proj_inplace_z (double *X, const int N)
{
    int n;
    _Complex double y;
    //struct timespec tic, toc;

    //Checks
    if (N<0) { fprintf(stderr,"error in proj_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);
    for (n=0; n<2*N; n+=2)
    {
        y = cproj(X[n]+1.0i*X[n+1]);
        memcpy(&X[n],(double *)&y,2*sizeof(double));
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
