//Gets square-root of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sqrt_s (float *Y, const float *X, const size_t N);
int sqrt_d (double *Y, const double *X, const size_t N);
int sqrt_c (float *Y, const float *X, const size_t N);
int sqrt_z (double *Y, const double *X, const size_t N);

int sqrt_inplace_s (float *X, const size_t N);
int sqrt_inplace_d (double *X, const size_t N);
int sqrt_inplace_c (float *X, const size_t N);
int sqrt_inplace_z (double *X, const size_t N);


int sqrt_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = sqrtf(*X); }

    return 0;
}


int sqrt_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = sqrt(*X); }
    
    return 0;
}


int sqrt_c (float *Y, const float *X, const size_t N)
{
    _Complex float y;

    for (size_t n=N; n>0u; --n, X+=2, ++Y)
    {
        y = csqrtf(*X + 1.0if**(X+1));
        *Y = *(float *)&y; *++Y = *((float *)&y+1);
    }
    
    return 0;
}


int sqrt_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=N; n>0u; --n, X+=2, ++Y)
    {
        y = csqrt(*X + 1.0i**(X+1));
        *Y = *(double *)&y; *++Y = *((double *)&y+1);
    }
    
    return 0;
}


int sqrt_inplace_s (float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = sqrtf(*X); }

    return 0;
}


int sqrt_inplace_d (double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = sqrt(*X); }
    
    return 0;
}


int sqrt_inplace_c (float *X, const size_t N)
{
    _Complex float y;

    for (size_t n=N; n>0u; --n, ++X)
    {
        y = csqrtf(*X + 1.0if**(X+1));
        *X = *(float *)&y; *++X = *((float *)&y+1);
    }
    
    return 0;
}


int sqrt_inplace_z (double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=N; n>0u; --n, ++X)
    {
        y = csqrt(*X + 1.0i**(X+1));
        *X = *(double *)&y; *++X = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
