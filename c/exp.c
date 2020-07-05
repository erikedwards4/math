//Gets exp of input X element-wise (e.^X).
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int exp_s (float *Y, const float *X, const size_t N);
int exp_d (double *Y, const double *X, const size_t N);
int exp_c (float *Y, const float *X, const size_t N);
int exp_z (double *Y, const double *X, const size_t N);

int exp_inplace_s (float *X, const size_t N);
int exp_inplace_d (double *X, const size_t N);
int exp_inplace_c (float *X, const size_t N);
int exp_inplace_z (double *X, const size_t N);


int exp_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = expf(X[n]); }

    return 0;
}


int exp_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = exp(X[n]); }
    
    return 0;
}


int exp_c (float *Y, const float *X, const size_t N)
{
    _Complex float y;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        y = cexpf(X[n2]+1.0if*X[n2+1]);
        memcpy(&Y[n2],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int exp_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        y = cexp(X[n2]+1.0i*X[n2+1]);
        memcpy(&Y[n2],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int exp_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = expf(X[n]); }

    return 0;
}


int exp_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = exp(X[n]); }
    
    return 0;
}


int exp_inplace_c (float *X, const size_t N)
{
    _Complex float x;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        x = cexpf(X[n2]+1.0if*X[n2+1]);
        memcpy(&X[n2],(float *)&x,2*sizeof(float));
    }
    
    return 0;
}


int exp_inplace_z (double *X, const size_t N)
{
    _Complex double x;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        x = cexp(X[n2]+1.0i*X[n2+1]);
        memcpy(&X[n2],(double *)&x,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
