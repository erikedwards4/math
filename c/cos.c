//Gets cosine (cos) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cos_s (float *Y, const float *X, const size_t N);
int cos_d (double *Y, const double *X, const size_t N);
int cos_c (float *Y, const float *X, const size_t N);
int cos_z (double *Y, const double *X, const size_t N);

int cos_inplace_s (float *X, const size_t N);
int cos_inplace_d (double *X, const size_t N);
int cos_inplace_c (float *X, const size_t N);
int cos_inplace_z (double *X, const size_t N);


int cos_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = cosf(*X); }

    return 0;
}


int cos_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = cos(*X); }
    
    return 0;
}


int cos_c (float *Y, const float *X, const size_t N)
{
    _Complex float y;
    
    for (size_t n=0; n<N; ++n, X+=2, ++Y)
    {
        y = ccosf(*X + 1.0if**(X+1));
        *Y = *(float *)&y; *++Y = *((float *)&y+1);
    }
    
    return 0;
}


int cos_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;
    
    for (size_t n=0; n<N; ++n, X+=2, ++Y)
    {
        y = ccos(*X + 1.0i**(X+1));
        *Y = *(double *)&y; *++Y = *((double *)&y+1);
    }
    
    return 0;
}


int cos_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = cosf(*X); }

    return 0;
}


int cos_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = cos(*X); }
    
    return 0;
}


int cos_inplace_c (float *X, const size_t N)
{
    _Complex float y;
    
    for (size_t n=0; n<N; ++n, ++X)
    {
        y = ccosf(*X + 1.0if**(X+1));
        *X = *(float *)&y; *++X = *((float *)&y+1);
    }
    
    return 0;
}


int cos_inplace_z (double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0; n<N; ++n, ++X)
    {
        y = ccos(*X + 1.0i**(X+1));
        *X = *(double *)&y; *++X = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
