//Gets log-gamma-function of each element of X.
//This has in-place and not-in-place versions.
//For complex input, gets lgamma of real and imag parts separately.

#include <stdio.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int lgamma_s (float *Y, const float *X, const size_t N);
int lgamma_d (double *Y, const double *X, const size_t N);
int lgamma_c (float *Y, const float *X, const size_t N);
int lgamma_z (double *Y, const double *X, const size_t N);

int lgamma_inplace_s (float *X, const size_t N);
int lgamma_inplace_d (double *X, const size_t N);
int lgamma_inplace_c (float *X, const size_t N);
int lgamma_inplace_z (double *X, const size_t N);


int lgamma_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = lgammaf(X[n]); }

    return 0;
}


int lgamma_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = lgamma(X[n]); }
    
    return 0;
}


int lgamma_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++) { Y[n] = lgammaf(X[n]); }
    
    return 0;
}


int lgamma_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++) { Y[n] = lgamma(X[n]); }
    
    return 0;
}


int lgamma_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = lgammaf(X[n]); }

    return 0;
}


int lgamma_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = lgamma(X[n]); }
    
    return 0;
}


int lgamma_inplace_c (float *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++) { X[n] = lgammaf(X[n]); }
    
    return 0;
}


int lgamma_inplace_z (double *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++) { X[n] = lgamma(X[n]); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
