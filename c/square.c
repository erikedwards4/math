//This just squares input X element-wise.
//For complex data, this is |X|.^2, i.e. Xr*Xr + Xi*Xi, and is sometimes called the (element-wise) norm.

//This has in-place and not-in-place versions.

#include <stdio.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int square_s (float *Y, const float *X, const size_t N);
int square_d (double *Y, const double *X, const size_t N);
int square_c (float *Y, const float *X, const size_t N);
int square_z (double *Y, const double *X, const size_t N);

int square_inplace_s (float *X, const size_t N);
int square_inplace_d (double *X, const size_t N);


int square_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X * *X; }

    return 0;
}


int square_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X * *X; }
    
    return 0;
}


int square_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, X+=2, ++Y) { *Y = *X**X + *(X+1)**(X+1); }
    
    return 0;
}


int square_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, X+=2, ++Y) { *Y = *X**X + *(X+1)**(X+1); }
    
    return 0;
}


int square_inplace_s (float *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X) { *X *= *X; }

    return 0;
}


int square_inplace_d (double *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X) { *X *= *X; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
