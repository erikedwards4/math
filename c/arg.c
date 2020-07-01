//Gets complex argument of input X element-wise.
//This is the phase angle in (-pi pi], equal to atan2(xi,xr).
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
//#include <complex.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int arg_s (float *Y, const float *X, const int N);
int arg_d (double *Y, const double *X, const int N);
int arg_c (float *Y, const float *X, const int N);
int arg_z (double *Y, const double *X, const int N);

int arg_inplace_s (float *X, const int N);
int arg_inplace_d (double *X, const int N);
int arg_inplace_c (float *X, const int N);
int arg_inplace_z (double *X, const int N);


int arg_s (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in arg_s: N (num elements X) must be nonnegative\n"); return 1; }

    //for (int n=0; n<N; n++) { Y[n] = atan2f(0.0f,X[n]); }
    for (int n=0; n<N; n++) { Y[n] = (X[n]<0.0f) ? (float)M_PI : 0.0f; }

    return 0;
}


int arg_d (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in arg_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { Y[n] = (X[n]<0.0) ? M_PI : 0.0; }
    
    return 0;
}


int arg_c (float *Y, const float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in arg_s: N (num elements X) must be nonnegative\n"); return 1; }

    int n = 0, n2 = 0;
    while (n<N) { Y[n] = atan2f(X[n2+1],X[n2]); n++; n2+=2; }
    //while (n<N) { Y[n] = cargf(X[n2]+1.0if*X[n2+1]); n++; n2+=2; }
    
    return 0;
}


int arg_z (double *Y, const double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in arg_z: N (num elements X) must be nonnegative\n"); return 1; }

    int n = 0, n2 = 0;
    while (n<N) { Y[n] = atan2(X[n2+1],X[n2]); n++; n2+=2; }
    
    return 0;
}


int arg_inplace_s (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in arg_inplace_s: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = (X[n]<0.0f) ? (float)M_PI : 0.0f; }
    //for (n=0; n<N; n++) { X[n] = (X[n]>0.0f)*(float)M_PI; }

    return 0;
}


int arg_inplace_d (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in arg_inplace_d: N (num elements X) must be nonnegative\n"); return 1; }

    for (int n=0; n<N; n++) { X[n] = (X[n]<0.0) ? M_PI : 0.0; }
    
    return 0;
}


int arg_inplace_c (float *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in arg_inplace_c: N (num elements X) must be nonnegative\n"); return 1; }

    int n = 0, n2 = 0;
    while (n<N) { X[n] = atan2f(X[n2+1],X[n2]); n++; n2+=2; }
    //while (n<N) { X[n] = cargf(X[n2]+1.0if*X[n2+1]); n++; n2+=2; }

    return 0;
}


int arg_inplace_z (double *X, const int N)
{
    if (N<0) { fprintf(stderr,"error in arg_inplace_z: N (num elements X) must be nonnegative\n"); return 1; }

    int n = 0, n2 = 0;
    while (n<N) { X[n] = atan2(X[n2+1],X[n2]); n++; n2+=2; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
