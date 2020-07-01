//This just squares input X element-wise.
//For complex data, this is |X|.^2, i.e. Xr*Xr + Xi*Xi.

//This has in-place and not-in-place versions.

#include <stdio.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int square_s (float *Y, const float *X, const int N);
int square_d (double *Y, const double *X, const int N);
int square_c (float *Y, const float *X, const int N);
int square_z (double *Y, const double *X, const int N);

int square_inplace_s (float *X, const int N);
int square_inplace_d (double *X, const int N);
int square_inplace_c (float *X, const int N);
int square_inplace_z (double *X, const int N);


int square_s (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in square_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = X[n]*X[n]; }

    return 0;
}


int square_d (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in square_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = X[n]*X[n]; }
    
    return 0;
}


int square_c (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in square_c: N (num elements X) must be nonnegative\n"); return 1; }

    int n = 0, n2 = 0;

    while (n<N) { Y[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1]; n++; n2+=2; }
    
    return 0;
}


int square_z (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in square_z: N (num elements X) must be nonnegative\n"); return 1; }

    int n = 0, n2 = 0;

    while (n<N) { Y[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1]; n++; n2+=2; }
    
    return 0;
}


int square_inplace_s (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in square_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] *= X[n]; }

    return 0;
}


int square_inplace_d (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in square_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] *= X[n]; }
    
    return 0;
}


int square_inplace_c (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in square_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    int n = 0, n2 = 0;
    
    for (n=0, n2=0; n<N; n++, n2+=2) { X[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1]; n++; n2+=2; }
    
    return 0;
}


int square_inplace_z (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in square_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    int n = 0, n2 = 0;

    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    //for (int n=0; n<N; n++) { X[n] = X[n*2]*X[n*2] + X[n*2+1]*X[n*2+1]; }
    while (n<N) { X[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1]; n++; n2+=2; }
    //for (int n=0; n<2*N; n++) { X[n] *= X[n]; }
    //for (int n=0; n<N; n++) { X[n] = X[2*n] + X[2*n+1]; }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
