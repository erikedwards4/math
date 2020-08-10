//Gets gamma-function (special function) of each element of X.
//This has in-place and not-in-place versions.

//For complex input, the ctgamma function is not usually available for complex.h.

#include <stdio.h>
#include <math.h>
#include <complex.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int tgamma_s (float *Y, const float *X, const size_t N);
int tgamma_d (double *Y, const double *X, const size_t N);
int tgamma_c (float *Y, const float *X, const size_t N);
int tgamma_z (double *Y, const double *X, const size_t N);

int tgamma_inplace_s (float *X, const size_t N);
int tgamma_inplace_d (double *X, const size_t N);
int tgamma_inplace_c (float *X, const size_t N);
int tgamma_inplace_z (double *X, const size_t N);


int tgamma_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = tgammaf(*X); }

    return 0;
}


int tgamma_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = tgamma(*X); }
    
    return 0;
}


// int tgamma_c (float *Y, const float *X, const size_t N)
// {
//     _Complex float y;

//     for (size_t n=0; n<N; ++n, X+=2, ++Y)
//     {
//         y = ctgammaf(*X + 1.0if**(X+1));
//         *Y = *(float *)&y; *++Y = *((float *)&y+1);
//     }
    
//     return 0;
// }


// int tgamma_z (double *Y, const double *X, const size_t N)
// {
//     _Complex double y;

//     for (size_t n=0; n<N; ++n, X+=2, ++Y)
//     {
//         y = ctgamma(*X + 1.0i**(X+1));
//         *Y = *(double *)&y; *++Y = *((double *)&y+1);
//     }
    
//     return 0;
// }


int tgamma_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = tgammaf(*X); }

    return 0;
}


int tgamma_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = tgamma(*X); }
    
    return 0;
}


// int tgamma_inplace_c (float *X, const size_t N)
// {
//     _Complex float y;

//     for (size_t n=0; n<N; ++n, ++X)
//     {
//         y = ctgamma(*X + 1.0if**(X+1));
//         *X = *(float *)&y; *++X = *((float *)&y+1);
//     }
    
//     return 0;
// }


// int tgamma_inplace_z (double *X, const size_t N)
// {
//     _Complex double y;

//     for (size_t n=0; n<N; ++n, ++X)
//     {
//         y = ctgamma(*X + 1.0i**(X+1));
//         *X = *(double *)&y; *++X = *((double *)&y+1);
//     }
    
//     return 0;
// }


#ifdef __cplusplus
}
}
#endif
