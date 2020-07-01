//Gets ceil of each element of X (rounds toward Inf).
//This has in-place and not-in-place versions.
//For complex input, ceils real and imag parts separately.

#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ceil_s (float *Y, const float *X, const int N);
int ceil_d (double *Y, const double *X, const int N);
int ceil_c (float *Y, const float *X, const int N);
int ceil_z (double *Y, const double *X, const int N);

int ceil_inplace_s (float *X, const int N);
int ceil_inplace_d (double *X, const int N);
int ceil_inplace_c (float *X, const int N);
int ceil_inplace_z (double *X, const int N);


int ceil_s (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in ceil_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = ceilf(X[n]); }

    return 0;
}


int ceil_d (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in ceil_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = ceil(X[n]); }
    
    return 0;
}


int ceil_c (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in ceil_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { Y[n] = ceilf(X[n]); }
    
    return 0;
}


int ceil_z (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in ceil_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { Y[n] = ceil(X[n]); }
    
    return 0;
}


int ceil_inplace_s (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in ceil_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = ceilf(X[n]); }

    return 0;
}


int ceil_inplace_d (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in ceil_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = ceil(X[n]); }
    
    return 0;
}


int ceil_inplace_c (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in ceil_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { X[n] = ceilf(X[n]); }
    
    return 0;
}


int ceil_inplace_z (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in ceil_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { X[n] = ceil(X[n]); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
