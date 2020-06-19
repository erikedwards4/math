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
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in sinh_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = sinhf(X[n]); }

    return 0;
}


int sinh_d (double *Y, const double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in sinh_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = sinh(X[n]); }
    
    return 0;
}


int sinh_c (float *Y, const float *X, const int N)
{
    int n2;
    _Complex float y;

    //Checks
    if (N<0) { fprintf(stderr,"error in sinh_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (n2=0; n2<2*N; n2+=2)
    {
        y = csinhf(X[n2]+1.0if*X[n2+1]);
        memcpy(&Y[n2],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int sinh_z (double *Y, const double *X, const int N)
{
    int n2;
    _Complex double y;

    //Checks
    if (N<0) { fprintf(stderr,"error in sinh_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (n2=0; n2<2*N; n2+=2)
    {
        y = csinh(X[n2]+1.0i*X[n2+1]);
        memcpy(&Y[n2],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int sinh_inplace_s (float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in sinh_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { X[n] = sinhf(X[n]); }

    return 0;
}


int sinh_inplace_d (double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in sinh_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { X[n] = sinh(X[n]); }
    
    return 0;
}


int sinh_inplace_c (float *X, const int N)
{
    int n2;
    _Complex float x;

    //Checks
    if (N<0) { fprintf(stderr,"error in sinh_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (n2=0; n2<2*N; n2+=2)
    {
        x = csinhf(X[n2]+1.0if*X[n2+1]);
        memcpy(&X[n2],(float *)&x,2*sizeof(float));
    }
    
    return 0;
}


int sinh_inplace_z (double *X, const int N)
{
    int n2;
    _Complex double x;

    //Checks
    if (N<0) { fprintf(stderr,"error in sinh_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (n2=0; n2<2*N; n2+=2)
    {
        x = csinh(X[n2]+1.0i*X[n2+1]);
        memcpy(&X[n2],(double *)&x,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
