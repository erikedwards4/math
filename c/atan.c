//Gets inverse tangent (atan) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int atan_s (float *Y, const float *X, const int N);
int atan_d (double *Y, const double *X, const int N);
int atan_c (float *Y, const float *X, const int N);
int atan_z (double *Y, const double *X, const int N);

int atan_inplace_s (float *X, const int N);
int atan_inplace_d (double *X, const int N);
int atan_inplace_c (float *X, const int N);
int atan_inplace_z (double *X, const int N);


int atan_s (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atan_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = atanf(X[n]); }

    return 0;
}


int atan_d (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atan_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = atan(X[n]); }
    
    return 0;
}


int atan_c (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atan_c: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex float y;

    for (int n2=0; n2<2*N; n2+=2)
    {
        y = catanf(X[n2]+1.0if*X[n2+1]);
        memcpy(&Y[n2],(float *)&y,2*sizeof(float));
    }
    
    return 0;
}


int atan_z (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atan_z: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex double y;

    for (int n2=0; n2<2*N; n2+=2)
    {
        y = catan(X[n2]+1.0i*X[n2+1]);
        memcpy(&Y[n2],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int atan_inplace_s (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atan_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = atanf(X[n]); }

    return 0;
}


int atan_inplace_d (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atan_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = atan(X[n]); }
    
    return 0;
}


int atan_inplace_c (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atan_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    _Complex float x;

    for (int n2=0; n2<2*N; n2+=2)
    {
        x = catanf(X[n2]+1.0if*X[n2+1]);
        memcpy(&X[n2],(float *)&x,2*sizeof(float));
    }
    
    return 0;
}


int atan_inplace_z (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in atan_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }
    
    _Complex double x;

    for (int n2=0; n2<2*N; n2+=2)
    {
        x = catan(X[n2]+1.0i*X[n2+1]);
        memcpy(&X[n2],(double *)&x,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
