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

int trunc_s (float *Y, const float *X, const size_t N);
int trunc_d (double *Y, const double *X, const size_t N);
int trunc_c (float *Y, const float *X, const size_t N);
int trunc_z (double *Y, const double *X, const size_t N);

int trunc_inplace_s (float *X, const size_t N);
int trunc_inplace_d (double *X, const size_t N);
int trunc_inplace_c (float *X, const size_t N);
int trunc_inplace_z (double *X, const size_t N);


int trunc_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = truncf(X[n]); }

    return 0;
}


int trunc_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = trunc(X[n]); }
    
    return 0;
}


int trunc_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++) { Y[n] = truncf(X[n]); }
    
    return 0;
}


int trunc_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++) { Y[n] = trunc(X[n]); }
    
    return 0;
}


int trunc_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = truncf(X[n]); }

    return 0;
}


int trunc_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = trunc(X[n]); }
    
    return 0;
}


int trunc_inplace_c (float *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++) { X[n] = truncf(X[n]); }
    
    return 0;
}


int trunc_inplace_z (double *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++) { X[n] = trunc(X[n]); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
