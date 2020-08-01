//Gets gamma-function (special function) of each element of X.
//This has in-place and not-in-place versions.
//For complex input, gets tgamma of real and imag parts separately.

#include <stdio.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int tgamma_s (float *Y, const float *X, const size_t N);
int tgamma_d (double *Y, const double *X, const size_t N);
int tgamma_c (float *Y, const float *X, const size_t N);
int tgamma_z (double *Y, const double *X, const size_t N);

int tgamma_inplace_s (float *X, const size_t N);
int tgamma_inplace_d (double *X, const size_t N);
int tgamma_inplace_c (float *X, const size_t N);
int tgamma_inplace_z (double *X, const size_t N);


int tgamma_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = tgammaf(*X); }

    return 0;
}


int tgamma_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = tgamma(*X); }
    
    return 0;
}


int tgamma_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = tgammaf(*X); }
    
    return 0;
}


int tgamma_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = tgamma(*X); }
    
    return 0;
}


int tgamma_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = tgammaf(*X); }

    return 0;
}


int tgamma_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = tgamma(*X); }
    
    return 0;
}


int tgamma_inplace_c (float *X, const size_t N)
{
    for (size_t n=0; n<2*N; ++n, ++X) { *X = tgammaf(*X); }
    
    return 0;
}


int tgamma_inplace_z (double *X, const size_t N)
{
    for (size_t n=0; n<2*N; ++n, ++X) { *X = tgamma(*X); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
