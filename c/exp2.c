//Gets exp2 of input X element-wise (2.^X).
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int exp2_s (float *Y, const float *X, const size_t N);
int exp2_d (double *Y, const double *X, const size_t N);
int exp2_c (float *Y, const float *X, const size_t N);
int exp2_z (double *Y, const double *X, const size_t N);

int exp2_inplace_s (float *X, const size_t N);
int exp2_inplace_d (double *X, const size_t N);
int exp2_inplace_c (float *X, const size_t N);
int exp2_inplace_z (double *X, const size_t N);


int exp2_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = exp2f(X[n]); }

    return 0;
}


int exp2_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = exp2(X[n]); }
    
    return 0;
}


int exp2_c (float *Y, const float *X, const size_t N)
{
    _Complex float y;

    for (size_t n=0; n<N; n++, X+=2)
    {
        y = cpowf(2.0f,*X+1.0if**(X+1));
        *Y++ = *(float *)&y; *Y++ = *((float *)&y+1);
    }
    
    return 0;
}


int exp2_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0; n<N; n++, X+=2)
    {
        y = cpow(2.0,*X+1.0i**(X+1));
        *Y++ = *(double *)&y; *Y++ = *((double *)&y+1);
    }
    
    return 0;
}


int exp2_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = exp2f(X[n]); }

    return 0;
}


int exp2_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = exp2(X[n]); }
    
    return 0;
}


int exp2_inplace_c (float *X, const size_t N)
{
    _Complex float y;

    for (size_t n=0; n<N; n++)
    {
        y = cpowf(2.0f,*X+1.0if**(X+1));
        *X++ = *(float *)&y; *X++ = *((float *)&y+1);
    }
    
    return 0;
}


int exp2_inplace_z (double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0; n<N; n++)
    {
        y = cpow(2.0,*X+1.0i**(X+1));
        *X++ = *(double *)&y; *X++ = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
