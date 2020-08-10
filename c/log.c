//Gets natural logarithm (ln) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int log_s (float *Y, const float *X, const size_t N);
int log_d (double *Y, const double *X, const size_t N);
int log_c (float *Y, const float *X, const size_t N);
int log_z (double *Y, const double *X, const size_t N);

int log_inplace_s (float *X, const size_t N);
int log_inplace_d (double *X, const size_t N);
int log_inplace_c (float *X, const size_t N);
int log_inplace_z (double *X, const size_t N);


int log_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = logf(*X); }

    return 0;
}


int log_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = log(*X); }
    
    return 0;
}


int log_c (float *Y, const float *X, const size_t N)
{
    _Complex float y;

    for (size_t n=0; n<N; ++n, X+=2, ++Y)
    {
        y = clogf(*X + 1.0if**(X+1));
        *Y = *(float *)&y; *++Y = *((float *)&y+1);
    }
    
    return 0;
}


int log_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0; n<N; ++n, X+=2, ++Y)
    {
        y = clog(*X + 1.0i**(X+1));
        *Y = *(double *)&y; *++Y = *((double *)&y+1);
    }
    
    return 0;
}


int log_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = logf(*X); }

    return 0;
}


int log_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = log(*X); }
    
    return 0;
}


int log_inplace_c (float *X, const size_t N)
{
    _Complex float y;

    for (size_t n=0; n<N; ++n, ++X)
    {
        y = clogf(*X + 1.0if**(X+1));
        *X = *(float *)&y; *++X = *((float *)&y+1);
    }
    
    return 0;
}


int log_inplace_z (double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0; n<N; ++n, ++X)
    {
        y = clog(*X + 1.0i**(X+1));
        *X = *(double *)&y; *++X = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
