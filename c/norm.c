//Gets norm (absolute value squared) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int norm_s (float *Y, const float *X, const size_t N);
int norm_d (double *Y, const double *X, const size_t N);
int norm_c (float *Y, const float *X, const size_t N);
int norm_z (double *Y, const double *X, const size_t N);

int norm_inplace_s (float *X, const size_t N);
int norm_inplace_d (double *X, const size_t N);
int norm_inplace_c (float *X, const size_t N);
int norm_inplace_z (double *X, const size_t N);


int norm_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = X[n]*X[n]; }

    return 0;
}


int norm_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = X[n]*X[n]; }
    
    return 0;
}


int norm_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++, X+=2) { *Y++ = *X**X + *(X+1)**(X+1); }
    
    return 0;
}


int norm_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<2*N; n++, X+=2) { *Y++ = *X**X + *(X+1)**(X+1); }
    
    return 0;
}


int norm_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = X[n]*X[n]; }

    return 0;
}


int norm_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = X[n]*X[n]; }
    
    return 0;
}


int norm_inplace_c (float *X, const size_t N)
{
    for (size_t n=0, n2=0; n<N; n++, n2+=2) { X[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1]; }
    
    return 0;
}


int norm_inplace_z (double *X, const size_t N)
{
    for (size_t n=0, n2=0; n<N; n++, n2+=2) { X[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1]; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
