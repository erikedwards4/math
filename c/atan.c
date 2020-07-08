//Gets inverse tangent (atan) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int atan_s (float *Y, const float *X, const size_t N);
int atan_d (double *Y, const double *X, const size_t N);
int atan_c (float *Y, const float *X, const size_t N);
int atan_z (double *Y, const double *X, const size_t N);

int atan_inplace_s (float *X, const size_t N);
int atan_inplace_d (double *X, const size_t N);
int atan_inplace_c (float *X, const size_t N);
int atan_inplace_z (double *X, const size_t N);


int atan_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = atanf(X[n]); }

    return 0;
}


int atan_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = atan(X[n]); }
    
    return 0;
}


int atan_c (float *Y, const float *X, const size_t N)
{
    _Complex float y;
    
    for (size_t n=0; n<N; n++, X+=2)
    {
        y = catanf(*X + 1.0if**(X+1));
        *Y++ = *(float *)&y; *Y++ = *((float *)&y+1);
    }
    
    return 0;
}


int atan_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;
    
    for (size_t n=0; n<N; n++, X+=2)
    {
        y = catan(*X + 1.0i**(X+1));
        *Y++ = *(double *)&y; *Y++ = *((double *)&y+1);
    }
    
    return 0;
}


int atan_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = atanf(X[n]); }

    return 0;
}


int atan_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = atan(X[n]); }
    
    return 0;
}


int atan_inplace_c (float *X, const size_t N)
{
    _Complex float y;
    
    for (size_t n=0; n<N; n++)
    {
        y = catanf(*X + 1.0if**(X+1));
        *X++ = *(float *)&y; *X++ = *((float *)&y+1);
    }
    
    return 0;
}


int atan_inplace_z (double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0; n<N; n++)
    {
        y = catan(*X + 1.0i**(X+1));
        *X++ = *(double *)&y; *X++ = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
