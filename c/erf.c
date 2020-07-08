//Gets error-function (special function) of each element of X.
//This has in-place and not-in-place versions.
//For complex input, gets erf of real and imag parts separately.

#include <stdio.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int erf_s (float *Y, const float *X, const size_t N);
int erf_d (double *Y, const double *X, const size_t N);
int erf_c (float *Y, const float *X, const size_t N);
int erf_z (double *Y, const double *X, const size_t N);

int erf_inplace_s (float *X, const size_t N);
int erf_inplace_d (double *X, const size_t N);
int erf_inplace_c (float *X, const size_t N);
int erf_inplace_z (double *X, const size_t N);


int erf_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = erff(X[n]); }

    return 0;
}


int erf_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = erf(X[n]); }
    
    return 0;
}


int erf_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++) { Y[n] = erff(X[n]); }
    
    return 0;
}


int erf_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++) { Y[n] = erf(X[n]); }
    
    return 0;
}


int erf_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = erff(X[n]); }

    return 0;
}


int erf_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = erf(X[n]); }
    
    return 0;
}


int erf_inplace_c (float *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++) { X[n] = erff(X[n]); }
    
    return 0;
}


int erf_inplace_z (double *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++) { X[n] = erf(X[n]); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
