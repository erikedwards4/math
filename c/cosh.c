//Gets hyberbolic cosine (cosh) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cosh_s (float *Y, const float *X, const size_t N);
int cosh_d (double *Y, const double *X, const size_t N);
int cosh_c (float *Y, const float *X, const size_t N);
int cosh_z (double *Y, const double *X, const size_t N);

int cosh_inplace_s (float *X, const size_t N);
int cosh_inplace_d (double *X, const size_t N);
int cosh_inplace_c (float *X, const size_t N);
int cosh_inplace_z (double *X, const size_t N);


int cosh_s (float *Y, const float *X, const size_t N)
{
    float xp, xm;

    //for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = coshf(*X); }
    for (size_t n=0; n<N; ++n, ++X, ++Y) { xp = expf(*X); xm = expf(-*X); *Y = 0.5f*(xp+xm); }

    return 0;
}


int cosh_d (double *Y, const double *X, const size_t N)
{
    double xp, xm;

    //for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = cosh(*X); }
    for (size_t n=0; n<N; ++n, ++X, ++Y) { xp = exp(*X); xm = exp(-*X); *Y = 0.5*(xp+xm); }
    
    return 0;
}


int cosh_c (float *Y, const float *X, const size_t N)
{
    _Complex float y, xp, xm;

    for (size_t n=0; n<N; ++n, X+=2, ++Y)
    {
        //y = ccoshf(X[n2]+1.0if*X[n2+1]);
        xp = cexpf(*X + 1.0if**(X+1));
        xm = cexpf(-*X - 1.0if**(X+1));
        y = 0.5f * (xp+xm);
        *Y++ = *(float *)&y; *Y = *((float *)&y+1);
    }
    
    return 0;
}


int cosh_z (double *Y, const double *X, const size_t N)
{
    _Complex double y, xp, xm;
    
    for (size_t n=0; n<N; ++n, X+=2, ++Y)
    {
        xp = cexp(*X + 1.0i**(X+1));
        xm = cexp(-*X - 1.0i**(X+1));
        y = 0.5 * (xp+xm);
        *Y++ = *(double *)&y; *Y = *((double *)&y+1);
    }
    
    return 0;
}


int cosh_inplace_s (float *X, const size_t N)
{
    float xp, xm;

    //for (size_t n=0; n<N; ++n, ++X) { *X = coshf(*X); }
    for (size_t n=0; n<N; ++n, ++X) { xp = expf(*X); xm = expf(-*X); *X = 0.5f*(xp+xm); }

    return 0;
}


int cosh_inplace_d (double *X, const size_t N)
{
    double xp, xm;

    //for (size_t n=0; n<N; ++n, ++X) { *X = cosh(*X); }
    for (size_t n=0; n<N; ++n, ++X) { xp = exp(*X); xm = exp(-*X); *X = 0.5*(xp+xm); }
    
    return 0;
}


int cosh_inplace_c (float *X, const size_t N)
{
    _Complex float y, xp, xm;

    for (size_t n=0; n<N; ++n, ++X)
    {
        //y = ccoshf(*X + 1.0if**(X+1));
        xp = cexpf(*X + 1.0if**(X+1));
        xm = cexpf(-*X - 1.0if**(X+1));
        y = 0.5f * (xp+xm);
        *X++ = *(float *)&y; *X = *((float *)&y+1);
    }
    
    return 0;
}


int cosh_inplace_z (double *X, const size_t N)
{
    _Complex double y, xp, xm;
    
    for (size_t n=0; n<N; ++n, ++X)
    {
        xp = cexp(*X + 1.0i**(X+1));
        xm = cexp(-*X - 1.0i**(X+1));
        y = 0.5 * (xp+xm);
        *X++ = *(double *)&y; *X = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
