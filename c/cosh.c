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

    //for (size_t n=0; n<N; n++) { Y[n] = coshf(X[n]); }
    for (size_t n=0; n<N; n++) { xp = expf(X[n]); xm = expf(-X[n]); Y[n] = 0.5f*(xp+xm); }

    return 0;
}


int cosh_d (double *Y, const double *X, const size_t N)
{
    double xp, xm;

    //for (size_t n=0; n<N; n++) { Y[n] = cosh(X[n]); }
    for (size_t n=0; n<N; n++) { xp = exp(X[n]); xm = exp(-X[n]); Y[n] = 0.5*(xp+xm); }
    
    return 0;
}


int cosh_c (float *Y, const float *X, const size_t N)
{
    _Complex float y, xp, xm;

    for (size_t n=0; n<N; n++, X+=2)
    {
        //y = ccoshf(X[n2]+1.0if*X[n2+1]);
        xp = cexpf(*X + 1.0if**(X+1));
        xm = cexpf(-*X - 1.0if**(X+1));
        y = 0.5f * (xp+xm);
        *Y++ = *(float *)&y; *Y++ = *((float *)&y+1);
    }
    
    return 0;
}


int cosh_z (double *Y, const double *X, const size_t N)
{
    _Complex double y, xp, xm;
    
    for (size_t n=0; n<N; n++, X+=2)
    {
        xp = cexp(*X + 1.0i**(X+1));
        xm = cexp(-*X - 1.0i**(X+1));
        y = 0.5 * (xp+xm);
        *Y++ = *(double *)&y; *Y++ = *((double *)&y+1);
    }
    
    return 0;
}


int cosh_inplace_s (float *X, const size_t N)
{
    float xp, xm;

    //for (size_t n=0; n<N; n++) { X[n] = coshf(X[n]); }
    for (size_t n=0; n<N; n++) { xp = expf(X[n]); xm = expf(-X[n]); X[n] = 0.5f*(xp+xm); }

    return 0;
}


int cosh_inplace_d (double *X, const size_t N)
{
    double xp, xm;

    //for (size_t n=0; n<N; n++) { X[n] = cosh(X[n]); }
    for (size_t n=0; n<N; n++) { xp = exp(X[n]); xm = exp(-X[n]); X[n] = 0.5*(xp+xm); }
    
    return 0;
}


int cosh_inplace_c (float *X, const size_t N)
{
    _Complex float y, xp, xm;

    for (size_t n=0; n<N; n++)
    {
        //y = ccoshf(*X + 1.0if**(X+1));
        xp = cexpf(*X + 1.0if**(X+1));
        xm = cexpf(-*X - 1.0if**(X+1));
        y = 0.5f * (xp+xm);
        *X++ = *(float *)&y; *X++ = *((float *)&y+1);
    }
    
    return 0;
}


int cosh_inplace_z (double *X, const size_t N)
{
    _Complex double y, xp, xm;
    
    for (size_t n=0; n<N; n++)
    {
        xp = cexp(*X + 1.0i**(X+1));
        xm = cexp(-*X - 1.0i**(X+1));
        y = 0.5 * (xp+xm);
        *X++ = *(double *)&y; *X++ = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
