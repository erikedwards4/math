//Gets log-gamma-function of each element of X.
//This has in-place and not-in-place versions.

//For complex input, the clgamma function is not usually available for complex.h.

#include <stdio.h>
#include <math.h>
#include <complex.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int lgamma_s (float *Y, const float *X, const size_t N);
int lgamma_d (double *Y, const double *X, const size_t N);
int lgamma_c (float *Y, const float *X, const size_t N);
int lgamma_z (double *Y, const double *X, const size_t N);

int lgamma_inplace_s (float *X, const size_t N);
int lgamma_inplace_d (double *X, const size_t N);
int lgamma_inplace_c (float *X, const size_t N);
int lgamma_inplace_z (double *X, const size_t N);


int lgamma_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = lgammaf(*X); }

    return 0;
}


int lgamma_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = lgamma(*X); }
    
    return 0;
}


// int lgamma_c (float *Y, const float *X, const size_t N)
// {
//     _Complex float y;

//     for (size_t n=0; n<N; ++n, X+=2, ++Y)
//     {
//         y = clgammaf(*X + 1.0if**(X+1));
//         *Y = *(float *)&y; *++Y = *((float *)&y+1);
//     }
    
//     return 0;
// }


// int lgamma_z (double *Y, const double *X, const size_t N)
// {
//     _Complex double y;

//     for (size_t n=0; n<N; ++n, X+=2, ++Y)
//     {
//         y = clgamma(*X + 1.0i**(X+1));
//         *Y = *(double *)&y; *++Y = *((double *)&y+1);
//     }
    
//     return 0;
// }


int lgamma_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = lgammaf(*X); }

    return 0;
}


int lgamma_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = lgamma(*X); }
    
    return 0;
}


// int lgamma_inplace_c (float *X, const size_t N)
// {
//     _Complex float y;

//     for (size_t n=0; n<N; ++n, ++X)
//     {
//         y = clgamma(*X + 1.0if**(X+1));
//         *X = *(float *)&y; *++X = *((float *)&y+1);
//     }
    
//     return 0;
// }


// int lgamma_inplace_z (double *X, const size_t N)
// {
//     _Complex double y;

//     for (size_t n=0; n<N; ++n, ++X)
//     {
//         y = clgamma(*X + 1.0i**(X+1));
//         *X = *(double *)&y; *++X = *((double *)&y+1);
//     }
    
//     return 0;
// }


#ifdef __cplusplus
}
}
#endif
