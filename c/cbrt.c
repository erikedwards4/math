//Gets cube-root of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cbrt_s (float *Y, const float *X, const size_t N);
int cbrt_d (double *Y, const double *X, const size_t N);
int cbrt_c (float *Y, const float *X, const size_t N);
int cbrt_z (double *Y, const double *X, const size_t N);

int cbrt_inplace_s (float *X, const size_t N);
int cbrt_inplace_d (double *X, const size_t N);
int cbrt_inplace_c (float *X, const size_t N);
int cbrt_inplace_z (double *X, const size_t N);


int cbrt_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = cbrtf(X[n]); }

    return 0;
}


int cbrt_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = cbrt(X[n]); }
    
    return 0;
}


int cbrt_c (float *Y, const float *X, const size_t N)
{
    const float p = 1.0f/3.0f;
    _Complex float y;

    for (size_t n=0; n<N; n++, X+=2)
    {
        y = cpowf(*X+1.0if**(X+1),p);
        *Y++ = *(float *)&y; *Y++ = *((float *)&y+1);
    }
    
    return 0;
}


int cbrt_z (double *Y, const double *X, const size_t N)
{
    const double p = 1.0/3.0;
    _Complex double y;

    for (size_t n=0; n<N; n++, X+=2)
    {
        y = cpow(*X+1.0i**(X+1),p);
        *Y++ = *(double *)&y; *Y++ = *((double *)&y+1);
    }
    
    return 0;
}


int cbrt_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = cbrtf(X[n]); }

    return 0;
}


int cbrt_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = cbrt(X[n]); }
    
    return 0;
}


int cbrt_inplace_c (float *X, const size_t N)
{
    const float p = 1.0f/3.0f;
    _Complex float y;

    for (size_t n=0; n<N; n++)
    {
        y = cpowf(*X+1.0if**(X+1),p);
        *X++ = *(float *)&y; *X++ = *((float *)&y+1);
    }
    
    return 0;
}


int cbrt_inplace_z (double *X, const size_t N)
{
    const double p = 1.0/3.0;
    _Complex double y;    

    for (size_t n=0; n<N; n++)
    {
        y = cpow(*X+1.0i**(X+1),p);
        *X++ = *(double *)&y; *X++ = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
