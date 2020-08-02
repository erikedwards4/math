//Gets inverse cosine (acos) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int acos_s (float *Y, const float *X, const size_t N);
int acos_d (double *Y, const double *X, const size_t N);
int acos_c (float *Y, const float *X, const size_t N);
int acos_z (double *Y, const double *X, const size_t N);

int acos_inplace_s (float *X, const size_t N);
int acos_inplace_d (double *X, const size_t N);
int acos_inplace_c (float *X, const size_t N);
int acos_inplace_z (double *X, const size_t N);


int acos_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = acosf(*X); }

    return 0;
}


int acos_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = acos(*X); }
    
    return 0;
}


int acos_c (float *Y, const float *X, const size_t N)
{
    _Complex float y;
    
    for (size_t n=0; n<N; ++n, X+=2, ++Y)
    {
        y = cacosf(*X + 1.0if**(X+1));
        *Y++ = *(float *)&y; *Y = *((float *)&y+1);
    }
    
    return 0;
}


int acos_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;
    
    for (size_t n=0; n<N; ++n, X+=2, ++Y)
    {
        y = cacos(*X + 1.0i**(X+1));
        *Y++ = *(double *)&y; *Y = *((double *)&y+1);
    }
    
    return 0;
}


int acos_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = acosf(*X); }

    return 0;
}


int acos_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = acos(*X); }
    
    return 0;
}


int acos_inplace_c (float *X, const size_t N)
{
    _Complex float y;
    
    for (size_t n=0; n<N; ++n, ++X)
    {
        y = cacosf(*X + 1.0if**(X+1));
        *X++ = *(float *)&y; *X = *((float *)&y+1);
    }
    
    return 0;
}


int acos_inplace_z (double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0; n<N; ++n, ++X)
    {
        y = cacos(*X + 1.0i**(X+1));
        *X++ = *(double *)&y; *X = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
