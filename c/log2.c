//Gets base-2 logarithm of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace openn {
extern "C" {
#endif

int log2_s (float *Y, const float *X, const int N);
int log2_d (double *Y, const double *X, const int N);
int log2_c (float *Y, const float *X, const int N);
int log2_z (double *Y, const double *X, const int N);

int log2_inplace_s (float *X, const int N);
int log2_inplace_d (double *X, const int N);
int log2_inplace_c (float *X, const int N);
int log2_inplace_z (double *X, const int N);


int log2_s (float *Y, const float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in log2_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = log2f(X[n]); }

    return 0;
}


int log2_d (double *Y, const double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in log2_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { Y[n] = log2(X[n]); }
    
    return 0;
}


int log2_c (float *Y, const float *X, const int N)
{
    int n2;
    _Complex float y;
    const _Complex float den = clogf(2.0f);

    //Checks
    if (N<0) { fprintf(stderr,"error in log2_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (n2=0; n2<2*N; n2+=2)
    {
        y = clogf(X[n2]+1.0if*X[n2+1]) / den;
        memcpy(&Y[n2],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int log2_z (double *Y, const double *X, const int N)
{
    int n2;
    _Complex double y;
    const _Complex double den = clog(2.0);

    //Checks
    if (N<0) { fprintf(stderr,"error in log2_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (n2=0; n2<2*N; n2+=2)
    {
        y = clog(X[n2]+1.0i*X[n2+1]) / den;
        memcpy(&Y[n2],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int log2_inplace_s (float *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in log2_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { X[n] = log2f(X[n]); }

    return 0;
}


int log2_inplace_d (double *X, const int N)
{
    int n;

    //Checks
    if (N<0) { fprintf(stderr,"error in log2_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (n=0; n<N; n++) { X[n] = log2(X[n]); }
    
    return 0;
}


int log2_inplace_c (float *X, const int N)
{
    int n2;
    _Complex float x;
    const _Complex float den = clogf(2.0f);

    //Checks
    if (N<0) { fprintf(stderr,"error in log2_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    for (n2=0; n2<2*N; n2+=2)
    {
        x = clogf(X[n2]+1.0if*X[n2+1]) / den;
        memcpy(&X[n2],(float *)&x,2*sizeof(float));
    }
    
    return 0;
}


int log2_inplace_z (double *X, const int N)
{
    int n2;
    _Complex double x;
    const _Complex double den = clog(2.0);

    //Checks
    if (N<0) { fprintf(stderr,"error in log2_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    for (n2=0; n2<2*N; n2+=2)
    {
        x = clog(X[n2]+1.0i*X[n2+1]) / den;
        memcpy(&X[n2],(double *)&x,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
