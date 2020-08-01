//Gets exp10 of input X element-wise (10.^X).
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int exp10_s (float *Y, const float *X, const size_t N);
int exp10_d (double *Y, const double *X, const size_t N);
int exp10_c (float *Y, const float *X, const size_t N);
int exp10_z (double *Y, const double *X, const size_t N);

int exp10_inplace_s (float *X, const size_t N);
int exp10_inplace_d (double *X, const size_t N);
int exp10_inplace_c (float *X, const size_t N);
int exp10_inplace_z (double *X, const size_t N);


int exp10_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = powf(10.0f,*X); }

    return 0;
}


int exp10_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = pow(10.0,*X); }
    
    return 0;
}


int exp10_c (float *Y, const float *X, const size_t N)
{
    _Complex float y;

    for (size_t n=0; n<N; ++n, X+=2)
    {
        y = cpowf(10.0f,*X+1.0if**(X+1));
        *Y++ = *(float *)&y; *Y++ = *((float *)&y+1);
    }
    
    return 0;
}


int exp10_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0; n<N; ++n, X+=2)
    {
        y = cpow(10.0,*X+1.0i**(X+1));
        *Y++ = *(double *)&y; *Y++ = *((double *)&y+1);
    }
    
    return 0;
}


int exp10_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = powf(10.0f,*X); }

    return 0;
}


int exp10_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = pow(10.0,*X); }
    
    return 0;
}


int exp10_inplace_c (float *X, const size_t N)
{
    _Complex float y;

    for (size_t n=0; n<N; ++n)
    {
        y = cpowf(10.0f,*X+1.0if**(X+1));
        *X++ = *(float *)&y; *X++ = *((float *)&y+1);
    }
    
    return 0;
}


int exp10_inplace_z (double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0; n<N; ++n)
    {
        y = cpow(10.0,*X+1.0i**(X+1));
        *X++ = *(double *)&y; *X++ = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
