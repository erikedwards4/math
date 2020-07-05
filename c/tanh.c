//Gets hyberbolic tangent (tanh) of input X element-wise.
// y = tanh(x) = (exp(x)-exp(-x))/(exp(x)+exp(-x))
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int tanh_s (float *Y, const float *X, const size_t N);
int tanh_d (double *Y, const double *X, const size_t N);
int tanh_c (float *Y, const float *X, const size_t N);
int tanh_z (double *Y, const double *X, const size_t N);

int tanh_inplace_s (float *X, const size_t N);
int tanh_inplace_d (double *X, const size_t N);
int tanh_inplace_c (float *X, const size_t N);
int tanh_inplace_z (double *X, const size_t N);


int tanh_s (float *Y, const float *X, const size_t N)
{
    float xp, xm;

    //for (size_t n=0; n<N; n++) { Y[n] = tanhf(X[n]); }
    for (size_t n=0; n<N; n++) { xp = expf(X[n]); xm = expf(-X[n]); Y[n] = (xp-xm)/(xp+xm); }

    return 0;
}


int tanh_d (double *Y, const double *X, const size_t N)
{
    double xp, xm;
    
    //for (size_t n=0; n<N; n++) { Y[n] = tanh(X[n]); }
    for (size_t n=0; n<N; n++) { xp = exp(X[n]); xm = exp(-X[n]); Y[n] = (xp-xm)/(xp+xm); }
    
    return 0;
}


int tanh_c (float *Y, const float *X, const size_t N)
{
    _Complex float y, xp, xm;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        //y = ctanhf(X[n2]+1.0if*X[n2+1]);
        xp = cexpf(X[n2]+1.0if*X[n2+1]);
        xm = cexpf(-X[n2]-1.0if*X[n2+1]);
        y = (xp-xm) / (xp+xm);
        memcpy(&Y[n2],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int tanh_z (double *Y, const double *X, const size_t N)
{
    _Complex double y, xp, xm;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        //y = ctanh(X[n2]+1.0i*X[n2+1]);
        xp = cexp(X[n2]+1.0i*X[n2+1]);
        xm = cexp(-X[n2]-1.0i*X[n2+1]);
        y = (xp-xm) / (xp+xm);
        memcpy(&Y[n2],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int tanh_inplace_s (float *X, const size_t N)
{
    float xp, xm;

    //for (size_t n=0; n<N; n++) { X[n] = tanhf(X[n]); }
    //for (size_t n=0; n<N; n++) { x = expf(2.0f*X[n]); X[n] = (x-1.0f)/(x+1.0f); }   //this is faster, but less numerical accuracy
    for (size_t n=0; n<N; n++) { xp = expf(X[n]); xm = expf(-X[n]); X[n] = (xp-xm)/(xp+xm); }

    return 0;
}


int tanh_inplace_d (double *X, const size_t N)
{
    double xp, xm;

    for (size_t n=0; n<N; n++) { X[n] = tanh(X[n]); }
    for (size_t n=0; n<N; n++) { xp = exp(X[n]); xm = exp(-X[n]); X[n] = (xp-xm)/(xp+xm); }
    
    return 0;
}


int tanh_inplace_c (float *X, const size_t N)
{
    _Complex float x, xp, xm;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        //x = ctanhf(X[n2]+1.0if*X[n2+1]);
        xp = cexpf(X[n2]+1.0if*X[n2+1]);
        xm = cexpf(-X[n2]-1.0if*X[n2+1]);
        x = (xp-xm) / (xp+xm);
        memcpy(&X[n2],(float *)&x,2*sizeof(float));
    }
    
    return 0;
}


int tanh_inplace_z (double *X, const size_t N)
{
    _Complex double x, xp, xm;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        //x = ctanh(X[n2]+1.0i*X[n2+1]);
        xp = cexp(X[n2]+1.0i*X[n2+1]);
        xm = cexp(-X[n2]-1.0i*X[n2+1]);
        x = (xp-xm) / (xp+xm);
        memcpy(&X[n2],(double *)&x,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
