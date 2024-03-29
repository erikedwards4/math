//Gets error-function (special function) of each element of X.
//This has in-place and not-in-place versions.

//For complex input, the cerf function is not usually available for complex.h,
//so I use libcerf from: https://jugit.fz-juelich.de/mlz/libcerf.
//However, this is rarely used, and messes up CFFI, so I comment out for now.

#include <stdio.h>
#include <math.h>
#include <complex.h>
//#include <cerf.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int erf_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = erff(*X); }

    return 0;
}


int erf_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = erf(*X); }
    
    return 0;
}


// int erf_c (float *Y, const float *X, const size_t N)
// {
//     _Complex double y;

//     for (size_t n=N; n>0u; --n, X+=2, ++Y)
//     {
//         y = cerf((double)*X + 1.0i*(double)*(X+1));
//         *Y = (float)*(double *)&y; *++Y = (float)*((double *)&y+1);
//     }
    
//     return 0;
// }


// int erf_z (double *Y, const double *X, const size_t N)
// {
//     _Complex double y;

//     for (size_t n=N; n>0u; --n, X+=2, ++Y)
//     {
//         y = cerf(*X + 1.0i**(X+1));
//         *Y = *(double *)&y; *++Y = *((double *)&y+1);
//     }
    
//     return 0;
// }


int erf_inplace_s (float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = erff(*X); }

    return 0;
}


int erf_inplace_d (double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X) { *X = erf(*X); }
    
    return 0;
}


// int erf_inplace_c (float *X, const size_t N)
// {
//     _Complex double y;

//     for (size_t n=N; n>0u; --n, ++X)
//     {
//         y = cerf((double)*X + 1.0i*(double)*(X+1));
//         *X = (float)*(double *)&y; *++X = (float)*((double *)&y+1);
//     }
    
//     return 0;
// }


// int erf_inplace_z (double *X, const size_t N)
// {
//     _Complex double y;

//     for (size_t n=N; n>0u; --n, ++X)
//     {
//         y = cerf(*X + 1.0i**(X+1));
//         *X = *(double *)&y; *++X = *((double *)&y+1);
//     }
    
//     return 0;
// }


#ifdef __cplusplus
}
}
#endif
