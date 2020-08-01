//Gets complementary error-function (special function) of each element of X.
//This has in-place and not-in-place versions.
//For complex input, gets erfc of real and imag parts separately.

#include <stdio.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int erfc_s (float *Y, const float *X, const size_t N);
int erfc_d (double *Y, const double *X, const size_t N);
int erfc_c (float *Y, const float *X, const size_t N);
int erfc_z (double *Y, const double *X, const size_t N);

int erfc_inplace_s (float *X, const size_t N);
int erfc_inplace_d (double *X, const size_t N);
int erfc_inplace_c (float *X, const size_t N);
int erfc_inplace_z (double *X, const size_t N);


int erfc_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = erfcf(*X); }

    return 0;
}


int erfc_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = erfc(*X); }
    
    return 0;
}


int erfc_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = erfcf(*X); }
    
    return 0;
}


int erfc_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = erfc(*X); }
    
    return 0;
}


int erfc_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = erfcf(*X); }

    return 0;
}


int erfc_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = erfc(*X); }
    
    return 0;
}


int erfc_inplace_c (float *X, const size_t N)
{
    for (size_t n=0; n<2*N; ++n, ++X) { *X = erfcf(*X); }
    
    return 0;
}


int erfc_inplace_z (double *X, const size_t N)
{
    for (size_t n=0; n<2*N; ++n, ++X) { *X = erfc(*X); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
