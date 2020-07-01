//Gets hyberbolic sine (sinh) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sinh_s (float *Y, const float *X, const int N);
int sinh_d (double *Y, const double *X, const int N);
int sinh_c (float *Y, const float *X, const int N);
int sinh_z (double *Y, const double *X, const int N);

int sinh_inplace_s (float *X, const int N);
int sinh_inplace_d (double *X, const int N);
int sinh_inplace_c (float *X, const int N);
int sinh_inplace_z (double *X, const int N);


int sinh_s (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in sinh_s: N (num elements X) must be nonnegative\n"); return 1; }

    float xp, xm;

    //for (int n=0; n<N; n++) { Y[n] = sinhf(X[n]); }
    for (int n=0; n<N; n++) { xp = expf(X[n]); xm = expf(-X[n]); Y[n] = 0.5f*(xp-xm); }

    return 0;
}


int sinh_d (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in sinh_d: N (num elements X) must be nonnegative\n"); return 1; }
    
    double xp, xm;

    //for (int n=0; n<N; n++) { Y[n] = sinh(X[n]); }
    for (int n=0; n<N; n++) { xp = exp(X[n]); xm = exp(-X[n]); Y[n] = 0.5*(xp-xm); }
    
    return 0;
}


int sinh_c (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in sinh_c: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex float y, xp, xm;

    for (int n2=0; n2<2*N; n2+=2)
    {
        //y = csinhf(X[n2]+1.0if*X[n2+1]);
        xp = cexpf(X[n2]+1.0if*X[n2+1]); xm = cexpf(-X[n2]-1.0if*X[n2+1]); y = 0.5f*(xp-xm);
        memcpy(&Y[n2],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int sinh_z (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in sinh_z: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex double y, xp, xm;

    for (int n2=0; n2<2*N; n2+=2)
    {
        //y = csinh(X[n2]+1.0i*X[n2+1]);
        xp = cexp(X[n2]+1.0i*X[n2+1]); xm = cexp(-X[n2]-1.0i*X[n2+1]); y = 0.5*(xp-xm);
        memcpy(&Y[n2],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int sinh_inplace_s (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in sinh_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    float xp, xm;

    //for (int n=0; n<N; n++) { X[n] = sinhf(X[n]); }
    for (int n=0; n<N; n++) { xp = expf(X[n]); xm = expf(-X[n]); X[n] = 0.5f*(xp-xm); }

    return 0;
}


int sinh_inplace_d (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in sinh_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    double xp, xm;

    //for (int n=0; n<N; n++) { X[n] = sinh(X[n]); }
    for (int n=0; n<N; n++) { xp = exp(X[n]); xm = exp(-X[n]); X[n] = 0.5*(xp-xm); }
    
    return 0;
}


int sinh_inplace_c (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in sinh_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex float x, xp, xm;

    for (int n2=0; n2<2*N; n2+=2)
    {
        //x = csinhf(X[n2]+1.0if*X[n2+1]);
        xp = cexpf(X[n2]+1.0if*X[n2+1]); xm = cexpf(-X[n2]-1.0if*X[n2+1]); x = 0.5f*(xp-xm);
        memcpy(&X[n2],(float *)&x,2*sizeof(float));
    }
    
    return 0;
}


int sinh_inplace_z (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in sinh_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }
    
    _Complex double x, xp, xm;

    for (int n2=0; n2<2*N; n2+=2)
    {
        //x = csinh(X[n2]+1.0i*X[n2+1]);
        xp = cexp(X[n2]+1.0i*X[n2+1]); xm = cexp(-X[n2]-1.0i*X[n2+1]); x = 0.5*(xp-xm);
        memcpy(&X[n2],(double *)&x,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
