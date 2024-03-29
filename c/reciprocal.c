//Gets reciprocal (1/x) of input X element-wise: Y = X.^-1.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int reciprocal_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = 1.0f / *X; }

    return 0;
}


int reciprocal_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = 1.0 / *X; }
    
    return 0;
}


int reciprocal_c (float *Y, const float *X, const size_t N)
{
    float xr, xi, sq;

    for (size_t n=N; n>0u; --n, ++X, ++Y)
    {
        xr = *X; xi = *++X;
        sq = xr*xr + xi*xi;
        *Y = xr/sq; *++Y = -xi/sq;
    }
    
    return 0;
}


int reciprocal_z (double *Y, const double *X, const size_t N)
{
    double xr, xi, sq;

    for (size_t n=N; n>0u; --n, ++X, ++Y)
    {
        xr = *X; xi = *++X;
        sq = xr*xr + xi*xi;
        *Y = xr/sq; *++Y = -xi/sq;
    }
    
    return 0;
}


int reciprocal_inplace_s (float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = 1.0f / *X; }

    return 0;
}


int reciprocal_inplace_d (double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = 1.0 / *X; }
    
    return 0;
}


int reciprocal_inplace_c (float *X, const size_t N)
{
    float xr, xi, sq;

    for (size_t n=N; n>0u; --n, ++X)
    {
        xr = *X; xi = *(X+1);
        sq = xr*xr + xi*xi;
        *X = xr/sq; *++X = -xi/sq;
    }
    
    return 0;
}


int reciprocal_inplace_z (double *X, const size_t N)
{
    double xr, xi, sq;

    for (size_t n=N; n>0u; --n, ++X)
    {
        xr = *X; xi = *(X+1);
        sq = xr*xr + xi*xi;
        *X = xr/sq; *++X = -xi/sq;
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
