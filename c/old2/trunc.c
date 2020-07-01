//Truncates of each element of X (removes decimal part, so rounds toward 0).
//This has in-place and not-in-place versions.
//For complex input, truncs real and imag parts separately.

#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int trunc_s (float *Y, const float *X, const int N);
int trunc_d (double *Y, const double *X, const int N);
int trunc_c (float *Y, const float *X, const int N);
int trunc_z (double *Y, const double *X, const int N);

int trunc_inplace_s (float *X, const int N);
int trunc_inplace_d (double *X, const int N);
int trunc_inplace_c (float *X, const int N);
int trunc_inplace_z (double *X, const int N);


int trunc_s (float *Y, const float *X, const int N)
{
    //Checks
    if (N<0) { fprintf(stderr,"error in trunc_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = truncf(X[n]); }

    return 0;
}


int trunc_d (double *Y, const double *X, const int N)
{
    //Checks
    if (N<0) { fprintf(stderr,"error in trunc_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = trunc(X[n]); }
    
    return 0;
}


int trunc_c (float *Y, const float *X, const int N)
{
    //Checks
    if (N<0) { fprintf(stderr,"error in trunc_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { Y[n] = truncf(X[n]); }
    
    return 0;
}


int trunc_z (double *Y, const double *X, const int N)
{
    //Checks
    if (N<0) { fprintf(stderr,"error in trunc_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { Y[n] = trunc(X[n]); }
    
    return 0;
}


int trunc_inplace_s (float *X, const int N)
{
    //Checks
    if (N<0) { fprintf(stderr,"error in trunc_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = truncf(X[n]); }

    return 0;
}


int trunc_inplace_d (double *X, const int N)
{
    //Checks
    if (N<0) { fprintf(stderr,"error in trunc_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = trunc(X[n]); }
    
    return 0;
}


int trunc_inplace_c (float *X, const int N)
{
    //Checks
    if (N<0) { fprintf(stderr,"error in trunc_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { X[n] = truncf(X[n]); }
    
    return 0;
}


int trunc_inplace_z (double *X, const int N)
{
    //Checks
    if (N<0) { fprintf(stderr,"error in trunc_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { X[n] = trunc(X[n]); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
