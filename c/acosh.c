//Gets inverse hyberbolic cosine (acosh) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int acosh_s (float *Y, const float *X, const size_t N);
int acosh_d (double *Y, const double *X, const size_t N);
int acosh_c (float *Y, const float *X, const size_t N);
int acosh_z (double *Y, const double *X, const size_t N);

int acosh_inplace_s (float *X, const size_t N);
int acosh_inplace_d (double *X, const size_t N);
int acosh_inplace_c (float *X, const size_t N);
int acosh_inplace_z (double *X, const size_t N);


int acosh_s (float *Y, const float *X, const size_t N)
{
    //for (size_t n=0; n<N; n++) { Y[n] = acoshf(X[n]); }
    for (size_t n=0; n<N; n++) { Y[n] = logf(X[n] + sqrtf(fmaf(X[n],X[n],-1.0f))); }

    return 0;
}


int acosh_d (double *Y, const double *X, const size_t N)
{
    //for (size_t n=0; n<N; n++) { Y[n] = acosh(X[n]); }
    for (size_t n=0; n<N; n++) { Y[n] = log(X[n] + sqrt(fma(X[n],X[n],-1.0))); }
    
    return 0;
}


int acosh_c (float *Y, const float *X, const size_t N)
{
    _Complex float x, y;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        //y = cacoshf(X[n2]+1.0if*X[n2+1]);
        x = X[n2] + 1.0if*X[n2+1];
        y = clogf(x + csqrtf(x*x-1.0f));
        memcpy(&Y[n2],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int acosh_z (double *Y, const double *X, const size_t N)
{
    _Complex double x, y;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        //y = cacosh(X[n2]+1.0i*X[n2+1]);
        x = X[n2] + 1.0i*X[n2+1];
        y = clog(x + csqrt(x*x-1.0));
        memcpy(&Y[n2],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int acosh_inplace_s (float *X, const size_t N)
{
    //for (size_t n=0; n<N; n++) { X[n] = acoshf(X[n]); }
    for (size_t n=0; n<N; n++) { X[n] = logf(X[n] + sqrtf(fmaf(X[n],X[n],-1.0f))); }

    return 0;
}


int acosh_inplace_d (double *X, const size_t N)
{
    //for (size_t n=0; n<N; n++) { X[n] = acosh(X[n]); }
    for (size_t n=0; n<N; n++) { X[n] = log(X[n] + sqrt(fma(X[n],X[n],-1.0))); }
    
    return 0;
}


int acosh_inplace_c (float *X, const size_t N)
{
    _Complex float x;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        //x = cacoshf(X[n2]+1.0if*X[n2+1]);
        x = X[n2] + 1.0if*X[n2+1];
        x = clogf(x + csqrtf(x*x-1.0f));
        memcpy(&X[n2],(float *)&x,2*sizeof(float));
    }
    
    return 0;
}


int acosh_inplace_z (double *X, const size_t N)
{
    _Complex double x;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        //x = cacosh(X[n2]+1.0i*X[n2+1]);
        x = X[n2] + 1.0i*X[n2+1];
        x = clog(x + csqrt(x*x-1.0));
        memcpy(&X[n2],(double *)&x,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
