//Gets complementary error-function (special function) of each element of X.
//This has in-place and not-in-place versions.

//For complex input, the cerfc function is not usually available for complex.h,
//so I use libcerf from: https://jugit.fz-juelich.de/mlz/libcerf.

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <cerf.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int erfc_s (float *Y, const float *X, const size_t N);
int erfc_d (double *Y, const double *X, const size_t N);
int erfc_c (float *Y, const float *X, const size_t N);
int erfc_z (double *Y, const double *X, const size_t N);

int erfc_inplace_s (float *X, const size_t N);
int erfc_inplace_d (double *X, const size_t N);
int erfc_inplace_c (float *X, const size_t N);
int erfc_inplace_z (double *X, const size_t N);


int erfc_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = erfcf(*X); }

    return 0;
}


int erfc_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = erfc(*X); }
    
    return 0;
}


int erfc_c (float *Y, const float *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0u; n<N; ++n, X+=2, ++Y)
    {
        y = cerfc((double)*X + 1.0i*(double)*(X+1));
        *Y = (float)*(double *)&y; *++Y = (float)*((double *)&y+1);
    }
    
    return 0;
}


int erfc_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0u; n<N; ++n, X+=2, ++Y)
    {
        y = cerfc((double)*X + 1.0i*(double)*(X+1));
        *Y = *(double *)&y; *++Y = *((double *)&y+1);
    }
    
    return 0;
}


int erfc_inplace_s (float *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X) { *X = erfcf(*X); }

    return 0;
}


int erfc_inplace_d (double *X, const size_t N)
{
    for (size_t n=0u; n<N; ++n, ++X) { *X = erfc(*X); }
    
    return 0;
}


int erfc_inplace_c (float *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0u; n<N; ++n, ++X)
    {
        y = cerfc((double)*X + 1.0i*(double)*(X+1));
        *X = (float)*(double *)&y; *++X = (float)*((double *)&y+1);
    }
    
    return 0;
}


int erfc_inplace_z (double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=0u; n<N; ++n, ++X)
    {
        y = cerfc(*X + 1.0i**(X+1));
        *X = *(double *)&y; *++X = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
