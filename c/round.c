//Rounds of each element of X (towards nearest integer).
//This has in-place and not-in-place versions.
//For complex input, rounds real and imag parts separately.

#include <stdio.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int round_s (float *Y, const float *X, const size_t N);
int round_d (double *Y, const double *X, const size_t N);
int round_c (float *Y, const float *X, const size_t N);
int round_z (double *Y, const double *X, const size_t N);

int round_inplace_s (float *X, const size_t N);
int round_inplace_d (double *X, const size_t N);
int round_inplace_c (float *X, const size_t N);
int round_inplace_z (double *X, const size_t N);


int round_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = roundf(*X); }

    return 0;
}


int round_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = round(*X); }
    
    return 0;
}


int round_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0u; n<2*N; ++n, ++X, ++Y) { *Y = roundf(*X); }
    
    return 0;
}


int round_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0u; n<2*N; ++n, ++X, ++Y) { *Y = round(*X); }
    
    return 0;
}


int round_inplace_s (float *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X) { *X = roundf(*X); }

    return 0;
}


int round_inplace_d (double *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X) { *X = round(*X); }
    
    return 0;
}


int round_inplace_c (float *X, const size_t N)
{
    for (size_t n=0u; n<2*N; ++n, ++X) { *X = roundf(*X); }
    
    return 0;
}


int round_inplace_z (double *X, const size_t N)
{
    for (size_t n=0u; n<2*N; ++n, ++X) { *X = round(*X); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
