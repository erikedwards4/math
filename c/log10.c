//Gets base-10 logarithm of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int log10_s (float *Y, const float *X, const size_t N);
int log10_d (double *Y, const double *X, const size_t N);
int log10_c (float *Y, const float *X, const size_t N);
int log10_z (double *Y, const double *X, const size_t N);

int log10_inplace_s (float *X, const size_t N);
int log10_inplace_d (double *X, const size_t N);
int log10_inplace_c (float *X, const size_t N);
int log10_inplace_z (double *X, const size_t N);


int log10_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = log10f(*X); }

    return 0;
}


int log10_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = log10(*X); }
    
    return 0;
}


int log10_c (float *Y, const float *X, const size_t N)
{
    const _Complex float den = clogf(10.0f);
    _Complex float y;

    for (size_t n=0u; n<N; ++n, X+=2, ++Y)
    {
        y = clogf(*X + 1.0if**(X+1)) / den;
        *Y = *(float *)&y; *++Y = *((float *)&y+1);
    }
    
    return 0;
}


int log10_z (double *Y, const double *X, const size_t N)
{
    const _Complex double den = clog(10.0);
    _Complex double y;

    for (size_t n=0u; n<N; ++n, X+=2, ++Y)
    {
        y = clog(*X + 1.0i**(X+1)) / den;
        *Y = *(double *)&y; *++Y = *((double *)&y+1);
    }
    
    return 0;
}


int log10_inplace_s (float *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X) { *X = log10f(*X); }

    return 0;
}


int log10_inplace_d (double *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X) { *X = log10(*X); }
    
    return 0;
}


int log10_inplace_c (float *X, const size_t N)
{
    const _Complex float den = clogf(10.0f);
    _Complex float y;

    for (size_t n=0u; n<N; ++n, ++X)
    {
        y = clogf(*X + 1.0if**(X+1)) / den;
        *X = *(float *)&y; *++X = *((float *)&y+1);
    }
    
    return 0;
}


int log10_inplace_z (double *X, const size_t N)
{
    const _Complex double den = clog(10.0);
    _Complex double y;

    for (size_t n=0u; n<N; ++n, ++X)
    {
        y = clog(*X + 1.0i**(X+1)) / den;
        *X = *(double *)&y; *++X = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
