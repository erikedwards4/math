//Gets inverse sine (asin) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int asin_s (float *Y, const float *X, const size_t N);
int asin_d (double *Y, const double *X, const size_t N);
int asin_c (float *Y, const float *X, const size_t N);
int asin_z (double *Y, const double *X, const size_t N);

int asin_inplace_s (float *X, const size_t N);
int asin_inplace_d (double *X, const size_t N);
int asin_inplace_c (float *X, const size_t N);
int asin_inplace_z (double *X, const size_t N);


int asin_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = asinf(*X); }

    return 0;
}


int asin_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = asin(*X); }
    
    return 0;
}


int asin_c (float *Y, const float *X, const size_t N)
{
    _Complex float y;
    
    for (size_t n=0; n<N; ++n, X+=2)
    {
        y = casinf(*X + 1.0if**(X+1));
        *Y++ = *(float *)&y; *Y++ = *((float *)&y+1);
    }
    
    return 0;
}


int asin_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;
    
    for (size_t n=0; n<N; ++n, X+=2)
    {
        y = casin(*X + 1.0i**(X+1));
        *Y++ = *(double *)&y; *Y++ = *((double *)&y+1);
    }
    
    return 0;
}


int asin_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = asinf(*X); }

    return 0;
}


int asin_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = asin(*X); }
    
    return 0;
}


int asin_inplace_c (float *X, const size_t N)
{
    _Complex float y;
    
    for (size_t n=0; n<N; ++n)
    {
        y = casinf(*X + 1.0if**(X+1));
        *X++ = *(float *)&y; *X++ = *((float *)&y+1);
    }
    
    return 0;
}


int asin_inplace_z (double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0; n<N; ++n)
    {
        y = casin(*X + 1.0i**(X+1));
        *X++ = *(double *)&y; *X++ = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
