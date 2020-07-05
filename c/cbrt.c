//Gets cube-root of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cbrt_s (float *Y, const float *X, const size_t N);
int cbrt_d (double *Y, const double *X, const size_t N);
int cbrt_c (float *Y, const float *X, const size_t N);
int cbrt_z (double *Y, const double *X, const size_t N);

int cbrt_inplace_s (float *X, const size_t N);
int cbrt_inplace_d (double *X, const size_t N);
int cbrt_inplace_c (float *X, const size_t N);
int cbrt_inplace_z (double *X, const size_t N);


int cbrt_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = cbrtf(X[n]); }

    return 0;
}


int cbrt_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = cbrt(X[n]); }
    
    return 0;
}


int cbrt_c (float *Y, const float *X, const size_t N)
{
    const float p = 1.0f/3.0f;
    _Complex float y;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        y = cpowf(X[n2]+1.0if*X[n2+1],p);
        memcpy(&Y[n2],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int cbrt_z (double *Y, const double *X, const size_t N)
{
    const double p = 1.0/3.0;
    _Complex double y;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        y = cpow(X[n2]+1.0i*X[n2+1],p);
        memcpy(&Y[n2],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int cbrt_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = cbrtf(X[n]); }

    return 0;
}


int cbrt_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = cbrt(X[n]); }
    
    return 0;
}


int cbrt_inplace_c (float *X, const size_t N)
{
    const float p = 1.0f/3.0f;
    _Complex float x;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        x = cpowf(X[n2]+1.0if*X[n2+1],p);
        memcpy(&X[n2],(float *)&x,2*sizeof(float));
    }
    
    return 0;
}


int cbrt_inplace_z (double *X, const size_t N)
{
    const double p = 1.0/3.0;
    _Complex double x;    

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        x = cpow(X[n2]+1.0i*X[n2+1],p);
        memcpy(&X[n2],(double *)&x,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
