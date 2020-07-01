//Gets floor of each element of X (rounds toward -Inf).
//This has in-place and not-in-place versions.
//For complex input, floors real and imag parts separately.

#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int floor_s (float *Y, const float *X, const int N);
int floor_d (double *Y, const double *X, const int N);
int floor_c (float *Y, const float *X, const int N);
int floor_z (double *Y, const double *X, const int N);

int floor_inplace_s (float *X, const int N);
int floor_inplace_d (double *X, const int N);
int floor_inplace_c (float *X, const int N);
int floor_inplace_z (double *X, const int N);


int floor_s (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in floor_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = floorf(X[n]); }

    return 0;
}


int floor_d (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in floor_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = floor(X[n]); }
    
    return 0;
}


int floor_c (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in floor_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { Y[n] = floorf(X[n]); }
    
    return 0;
}


int floor_z (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in floor_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { Y[n] = floor(X[n]); }
    
    return 0;
}


int floor_inplace_s (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in floor_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = floorf(X[n]); }

    return 0;
}


int floor_inplace_d (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in floor_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = floor(X[n]); }
    
    return 0;
}


int floor_inplace_c (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in floor_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { X[n] = floorf(X[n]); }
    
    return 0;
}


int floor_inplace_z (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in floor_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { X[n] = floor(X[n]); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
