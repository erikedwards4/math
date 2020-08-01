//Gets complex argument of input X element-wise.
//This is the phase angle in (-pi pi], equal to atan2(xi,xr).
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
//#include <complex.h>
//#include <time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int arg_s (float *Y, const float *X, const size_t N);
int arg_d (double *Y, const double *X, const size_t N);
int arg_c (float *Y, const float *X, const size_t N);
int arg_z (double *Y, const double *X, const size_t N);

int arg_inplace_s (float *X, const size_t N);
int arg_inplace_d (double *X, const size_t N);
int arg_inplace_c (float *X, const size_t N);
int arg_inplace_z (double *X, const size_t N);


int arg_s (float *Y, const float *X, const size_t N)
{
    //for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = atan2f(0.0f,*X); }
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = (float)(*X<0.0f) * (float)M_PI; }

    return 0;
}


int arg_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = (double)(*X<0.0) * M_PI; }
    
    return 0;
}


int arg_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, X+=2, ++Y) { *Y = atan2f(*(X+1),*X); }
    //for (size_t n=0, n2=0; n<N; ++n, n2+=2) { *Y = cargf(X[n2]+1.0if*X[n2+1]); }
    
    return 0;
}


int arg_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, X+=2, ++Y) { *Y = atan2(*(X+1),*X); }
    
    return 0;
}


int arg_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = (float)(*X<0.0f) * (float)M_PI; }
    //for (size_t n=0; n<N; ++n) { *X = (*X>0.0f)*(float)M_PI; }

    return 0;
}


int arg_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = (double)(*X<0.0) * M_PI; }
    
    return 0;
}


int arg_inplace_c (float *X, const size_t N)
{
    for (size_t n=0, n2=0; n<N; ++n, n2+=2) { X[n] = atan2f(X[n2+1],X[n2]); }
    //for (size_t n=0, n2=0; n<N; ++n, n2+=2) { X[n] = cargf(X[n2]+1.0if*X[n2+1]); }

    return 0;
}


int arg_inplace_z (double *X, const size_t N)
{
    for (size_t n=0, n2=0; n<N; ++n, n2+=2) { X[n] = atan2(X[n2+1],X[n2]); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
