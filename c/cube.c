//This just cubes input X element-wise.
//For complex data, this is |X|.^2, i.e. Xr*Xr + Xi*Xi.

//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cube_s (float *Y, const float *X, const int N);
int cube_d (double *Y, const double *X, const int N);
int cube_c (float *Y, const float *X, const int N);
int cube_z (double *Y, const double *X, const int N);

int cube_inplace_s (float *X, const int N);
int cube_inplace_d (double *X, const int N);
int cube_inplace_c (float *X, const int N);
int cube_inplace_z (double *X, const int N);


int cube_s (float *Y, const float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in cube_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = X[n]*X[n]*X[n]; }

    return 0;
}


int cube_d (double *Y, const double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in cube_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = X[n]*X[n]*X[n]; }
    
    return 0;
}


int cube_c (float *Y, const float *X, const int N)
{
    int n, n2;

    //Checks
    if (N<0) { fprintf(stderr,"error in cube_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0, n2=0; n<N; n++, n2+=2)
    {
        Y[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1];
        Y[n] *= sqrtf(Y[n]);
    }
    
    return 0;
}


int cube_z (double *Y, const double *X, const int N)
{
    int n, n2;

    //Checks
    if (N<0) { fprintf(stderr,"error in cube_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0, n2=0; n<N; n++, n2+=2)
    {
        Y[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1];
        Y[n] *= sqrt(Y[n]);
    }
    
    return 0;
}


int cube_inplace_s (float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in cube_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { X[n] = X[n]*X[n]*X[n]; }

    return 0;
}


int cube_inplace_d (double *X, const int N)
{
    int n;
    //struct timespec tic, toc;

    //Checks
    if (N<0) { fprintf(stderr,"error in cube_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);
    for (n=0; n<N; n++) { X[n] = X[n]*X[n]*X[n]; }
    //for (n=0; n<N; n++) { X[n] *= X[n]*X[n]; }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int cube_inplace_c (float *X, const int N)
{
    int n, n2;

    //Checks
    if (N<0) { fprintf(stderr,"error in cube_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0, n2=0; n<N; n++, n2+=2)
    {
        X[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1];
        X[n] *= sqrtf(X[n]);
    }
    
    return 0;
}


int cube_inplace_z (double *X, const int N)
{
    int n, n2;
    //struct timespec tic, toc;

    //Checks
    if (N<0) { fprintf(stderr,"error in cube_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&tic);
    for (n=0, n2=0; n<N; n++, n2+=2)
    {
        X[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1];
        X[n] *= sqrt(X[n]);
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
