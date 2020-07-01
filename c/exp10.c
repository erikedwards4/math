//Gets exp10 of input X element-wise (10.^X).
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int exp10_s (float *Y, const float *X, const int N);
int exp10_d (double *Y, const double *X, const int N);
int exp10_c (float *Y, const float *X, const int N);
int exp10_z (double *Y, const double *X, const int N);

int exp10_inplace_s (float *X, const int N);
int exp10_inplace_d (double *X, const int N);
int exp10_inplace_c (float *X, const int N);
int exp10_inplace_z (double *X, const int N);


int exp10_s (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in exp10_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = powf(10.0f,X[n]); }

    return 0;
}


int exp10_d (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in exp10_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = pow(10.0,X[n]); }
    
    return 0;
}


int exp10_c (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in exp10_c: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex float y;

    for (int n2=0; n2<2*N; n2+=2)
    {
        y = cpowf(10.0f,X[n2]+1.0if*X[n2+1]);
        memcpy(&Y[n2],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int exp10_z (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in exp10_z: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex double y;

    for (int n2=0; n2<2*N; n2+=2)
    {
        y = cpow(10.0,X[n2]+1.0i*X[n2+1]);
        memcpy(&Y[n2],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int exp10_inplace_s (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in exp10_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = powf(10.0f,X[n]); }

    return 0;
}


int exp10_inplace_d (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in exp10_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = pow(10.0,X[n]); }
    
    return 0;
}


int exp10_inplace_c (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in exp10_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex float x;

    for (int n2=0; n2<2*N; n2+=2)
    {
        x = cpowf(10.0f,X[n2]+1.0if*X[n2+1]);
        memcpy(&X[n2],(float *)&x,2*sizeof(float));
    }
    
    return 0;
}


int exp10_inplace_z (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in exp10_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }
    
    _Complex double x;    

    for (int n2=0; n2<2*N; n2+=2)
    {
        x = cpow(10.0,X[n2]+1.0i*X[n2+1]);
        memcpy(&X[n2],(double *)&x,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
