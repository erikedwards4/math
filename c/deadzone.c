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
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = (float)((*X>delta)-(*X<-delta)); }

    return 0;
}


int deadzone_d (double *Y, const double *X, const size_t N, const double delta)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = (double)((*X>delta)-(*X<-delta)); }
    
    return 0;
}


int deadzone_inplace_s (float *X, const size_t N, const float delta)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = (float)((*X>delta)-(*X<-delta)); }

    return 0;
}


int deadzone_inplace_d (double *X, const size_t N, const double delta)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = (double)((*X>delta)-(*X<-delta)); }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
