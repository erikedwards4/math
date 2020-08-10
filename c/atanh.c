//Gets inverse hyberbolic tangent (atanh) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int atanh_s (float *Y, const float *X, const size_t N);
int atanh_d (double *Y, const double *X, const size_t N);
int atanh_c (float *Y, const float *X, const size_t N);
int atanh_z (double *Y, const double *X, const size_t N);

int atanh_inplace_s (float *X, const size_t N);
int atanh_inplace_d (double *X, const size_t N);
int atanh_inplace_c (float *X, const size_t N);
int atanh_inplace_z (double *X, const size_t N);


int atanh_s (float *Y, const float *X, const size_t N)
{
    //for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = atanhf(*X); }
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = 0.5f * logf((1.0f+*X)/(1.0f-*X)); }

    return 0;
}


int atanh_d (double *Y, const double *X, const size_t N)
{
    //for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = atanh(*X); }
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = 0.5 * log((1.0+*X)/(1.0-*X)); }
    
    return 0;
}


int atanh_c (float *Y, const float *X, const size_t N)
{
    _Complex float x, y;

    for (size_t n=0; n<N; ++n, X+=2, ++Y)
    {
        //y = catanhf(X[n2]+1.0if*X[n2+1]);
        x = *X + 1.0if**(X+1);
        y = 0.5f + clogf((1.0f+x)/(1.0f-x));
        *Y = *(float *)&y; *++Y = *((float *)&y+1);
    }
    
    return 0;
}


int atanh_z (double *Y, const double *X, const size_t N)
{
    _Complex double x, y;

    for (size_t n=0; n<N; ++n, X+=2, ++Y)
    {
        //y = catanh(X[n2]+1.0if*X[n2+1]);
        x = *X + 1.0i**(X+1);
        y = 0.5 + clog((1.0+x)/(1.0-x));
        *Y = *(double *)&y; *++Y = *((double *)&y+1);
    }
    
    return 0;
}


int atanh_inplace_s (float *X, const size_t N)
{
    //for (size_t n=0; n<N; ++n, ++X) { *X = atanhf(*X); }
    for (size_t n=0; n<N; ++n, ++X) { *X = 0.5f * logf((1.0f+*X)/(1.0f-*X)); }

    return 0;
}


int atanh_inplace_d (double *X, const size_t N)
{
    //for (size_t n=0; n<N; ++n, ++X) { *X = atanh(*X); }
    for (size_t n=0; n<N; ++n, ++X) { *X = 0.5 * log((1.0+*X)/(1.0-*X)); }
    
    return 0;
}


int atanh_inplace_c (float *X, const size_t N)
{
    _Complex float x, y;

    for (size_t n=0; n<N; ++n, ++X)
    {
        //y = catanhf(X[n2]+1.0if*X[n2+1]);
        x = *X + 1.0if**(X+1);
        y = 0.5f + clogf((1.0f+x)/(1.0f-x));
        *X = *(float *)&y; *++X = *((float *)&y+1);
    }
    
    return 0;
}


int atanh_inplace_z (double *X, const size_t N)
{
    _Complex double x, y;

    for (size_t n=0; n<N; ++n, ++X)
    {
        //y = catanh(X[n2]+1.0i*X[n2+1]);
        x = *X + 1.0i**(X+1);
        y = 0.5 + clog((1.0+x)/(1.0-x));
        *X = *(double *)&y; *++X = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
