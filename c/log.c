//Gets natural logarithm (ln) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>
//#include <time.h>

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
    for (size_t n=0; n<N; n++) { Y[n] = logf(X[n]); }

    return 0;
}


int log_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = log(X[n]); }
    
    return 0;
}


int log_c (float *Y, const float *X, const size_t N)
{
    _Complex float y;

    for (size_t n=0; n<N; n++, X+=2)
    {
        y = clogf(*X + 1.0if**(X+1));
        *Y++ = *(float *)&y; *Y++ = *((float *)&y+1);
    }
    
    return 0;
}


int log_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0; n<N; n++, X+=2)
    {
        y = clog(*X + 1.0i**(X+1));
        *Y++ = *(double *)&y; *Y++ = *((double *)&y+1);
    }
    
    return 0;
}


int log_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = logf(X[n]); }

    return 0;
}


int log_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = log(X[n]); }
    
    return 0;
}


int log_inplace_c (float *X, const size_t N)
{
    _Complex float y;

    for (size_t n=0; n<N; n++)
    {
        y = clogf(*X + 1.0if**(X+1));
        *X++ = *(float *)&y; *X++ = *((float *)&y+1);
    }
    
    return 0;
}


int log_inplace_z (double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0; n<N; n++)
    {
        y = clog(*X + 1.0i**(X+1));
        *X++ = *(double *)&y; *X++ = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
