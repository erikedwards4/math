//Rounds of each element of X (towards nearest integer).
//This has in-place and not-in-place versions.
//For complex input, rounds real and imag parts separately.

#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int round_s (float *Y, const float *X, const int N);
int round_d (double *Y, const double *X, const int N);
int round_c (float *Y, const float *X, const int N);
int round_z (double *Y, const double *X, const int N);

int round_inplace_s (float *X, const int N);
int round_inplace_d (double *X, const int N);
int round_inplace_c (float *X, const int N);
int round_inplace_z (double *X, const int N);


int round_s (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in round_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = roundf(X[n]); }

    return 0;
}


int round_d (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in round_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = round(X[n]); }
    
    return 0;
}


int round_c (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in round_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { Y[n] = roundf(X[n]); }
    
    return 0;
}


int round_z (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in round_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { Y[n] = round(X[n]); }
    
    return 0;
}


int round_inplace_s (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in round_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = roundf(X[n]); }

    return 0;
}


int round_inplace_d (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in round_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = round(X[n]); }
    
    return 0;
}


int round_inplace_c (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in round_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { X[n] = roundf(X[n]); }
    
    return 0;
}


int round_inplace_z (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in round_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<2*N; n++) { X[n] = round(X[n]); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
