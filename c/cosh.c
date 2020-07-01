//Gets hyberbolic cosine (cosh) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cosh_s (float *Y, const float *X, const int N);
int cosh_d (double *Y, const double *X, const int N);
int cosh_c (float *Y, const float *X, const int N);
int cosh_z (double *Y, const double *X, const int N);

int cosh_inplace_s (float *X, const int N);
int cosh_inplace_d (double *X, const int N);
int cosh_inplace_c (float *X, const int N);
int cosh_inplace_z (double *X, const int N);


int cosh_s (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in cosh_s: N (num elements X) must be nonnegative\n"); return 1; }

    float xp, xm;

    //for (int n=0; n<N; n++) { Y[n] = coshf(X[n]); }
    for (int n=0; n<N; n++) { xp = expf(X[n]); xm = expf(-X[n]); Y[n] = 0.5f*(xp+xm); }

    return 0;
}


int cosh_d (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in cosh_d: N (num elements X) must be nonnegative\n"); return 1; }

    double xp, xm;

    //for (int n=0; n<N; n++) { Y[n] = cosh(X[n]); }
    for (int n=0; n<N; n++) { xp = exp(X[n]); xm = exp(-X[n]); Y[n] = 0.5*(xp+xm); }
    
    return 0;
}


int cosh_c (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in cosh_c: N (num elements X) must be nonnegative\n"); return 1; }
    
    _Complex float y, xp, xm;

    for (int n2=0; n2<2*N; n2+=2)
    {
        //y = ccoshf(X[n2]+1.0if*X[n2+1]);
        xp = cexpf(X[n2]+1.0if*X[n2+1]); xm = cexpf(-X[n2]-1.0if*X[n2+1]); y = 0.5f*(xp+xm);
        memcpy(&Y[n2],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int cosh_z (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in cosh_z: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex double y, xp, xm;

    for (int n2=0; n2<2*N; n2+=2)
    {
        //y = ccosh(X[n2]+1.0i*X[n2+1]);
        xp = cexp(X[n2]+1.0i*X[n2+1]); xm = cexp(-X[n2]-1.0i*X[n2+1]); y = 0.5*(xp+xm);
        memcpy(&Y[n2],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int cosh_inplace_s (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in cosh_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    float xp, xm;

    //for (int n=0; n<N; n++) { X[n] = coshf(X[n]); }
    for (int n=0; n<N; n++) { xp = expf(X[n]); xm = expf(-X[n]); X[n] = 0.5f*(xp+xm); }

    return 0;
}


int cosh_inplace_d (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in cosh_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    double xp, xm;

    //for (int n=0; n<N; n++) { X[n] = cosh(X[n]); }
    for (int n=0; n<N; n++) { xp = exp(X[n]); xm = exp(-X[n]); X[n] = 0.5*(xp+xm); }
    
    return 0;
}


int cosh_inplace_c (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in cosh_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex float x, xp, xm;

    for (int n2=0; n2<2*N; n2+=2)
    {
        //x = ccoshf(X[n2]+1.0if*X[n2+1]);
        xp = cexpf(X[n2]+1.0if*X[n2+1]); xm = cexpf(-X[n2]-1.0if*X[n2+1]); x = 0.5f*(xp+xm);
        memcpy(&X[n2],(float *)&x,2*sizeof(float));
    }
    
    return 0;
}


int cosh_inplace_z (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in cosh_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }
    
    _Complex double x, xp, xm;

    for (int n2=0; n2<2*N; n2+=2)
    {
        //x = ccosh(X[n2]+1.0i*X[n2+1]);
        xp = cexp(X[n2]+1.0i*X[n2+1]); xm = cexp(-X[n2]-1.0i*X[n2+1]); x = 0.5*(xp+xm);
        memcpy(&X[n2],(double *)&x,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
