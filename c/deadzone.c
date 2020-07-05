//This gets hard limiter with deadzone for each element of X.
//The width of the deadzone is from -delta to delta.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int deadzone_s (float *Y, const float *X, const size_t N, const float delta);
int deadzone_d (double *Y, const double *X, const size_t N, const double delta);

int deadzone_inplace_s (float *X, const size_t N, const float delta);
int deadzone_inplace_d (double *X, const size_t N, const double delta);


int deadzone_s (float *Y, const float *X, const size_t N, const float delta)
{
    for (size_t n=0; n<N; n++) { Y[n] = (X[n]>delta) - (X[n]<-delta); }

    return 0;
}


int deadzone_d (double *Y, const double *X, const size_t N, const double delta)
{
    for (size_t n=0; n<N; n++) { Y[n] = (X[n]>delta) - (X[n]<-delta); }
    
    return 0;
}


int deadzone_inplace_s (float *X, const size_t N, const float delta)
{
    for (size_t n=0; n<N; n++) { X[n] = (X[n]>delta) - (X[n]<-delta); }

    return 0;
}


int deadzone_inplace_d (double *X, const size_t N, const double delta)
{
    for (size_t n=0; n<N; n++) { X[n] = (X[n]>delta) - (X[n]<-delta); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
