//Gets hyberbolic sine (sinh) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int sinh_s (float *Y, const float *X, const size_t N)
{
    float xp, xm;

    //for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = sinhf(*X); }
    for (size_t n=N; n>0u; --n, ++X, ++Y) { xp = expf(*X); xm = expf(-*X); *Y = 0.5f*(xp-xm); }

    return 0;
}


int sinh_d (double *Y, const double *X, const size_t N)
{
    double xp, xm;

    //for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = sinh(*X); }
    for (size_t n=N; n>0u; --n, ++X, ++Y) { xp = exp(*X); xm = exp(-*X); *Y = 0.5*(xp-xm); }
    
    return 0;
}


int sinh_c (float *Y, const float *X, const size_t N)
{
    _Complex float y, xp, xm;
    
    for (size_t n=N; n>0u; --n, X+=2, ++Y)
    {
        //y = csinhf(*X + 1.0if**(X+1));
        xp = cexpf(*X + 1.0if**(X+1));
        xm = cexpf(-*X - 1.0if**(X+1));
        y = 0.5f * (xp-xm);
        //y = *X + 1.0if**(X+1);
        //y = 0.5f * (cexpf(y)-cexpf(-y));
        *Y = *(float *)&y; *++Y = *((float *)&y+1);
    }
    
    return 0;
}


int sinh_z (double *Y, const double *X, const size_t N)
{
    _Complex double y, xp, xm;
    
    for (size_t n=N; n>0u; --n, X+=2, ++Y)
    {
        xp = cexp(*X + 1.0i**(X+1));
        xm = cexp(-*X - 1.0i**(X+1));
        y = 0.5 * (xp-xm);
        *Y = *(double *)&y; *++Y = *((double *)&y+1);
    }
    
    return 0;
}


int sinh_inplace_s (float *X, const size_t N)
{
    float xp, xm;

    //for (size_t n=N; n>0u; --n, ++X) { *X = sinhf(*X); }
    for (size_t n=N; n>0u; --n, ++X) { xp = expf(*X); xm = expf(-*X); *X = 0.5f*(xp-xm); }

    return 0;
}


int sinh_inplace_d (double *X, const size_t N)
{
    double xp, xm;

    //for (size_t n=N; n>0u; --n, ++X) { *X = sinh(*X); }
    for (size_t n=N; n>0u; --n, ++X) { xp = exp(*X); xm = exp(-*X); *X = 0.5*(xp-xm); }
    
    return 0;
}


int sinh_inplace_c (float *X, const size_t N)
{
    _Complex float y, xp, xm;
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    for (size_t n=N; n>0u; --n, ++X)
    {
        //y = csinhf(*X + 1.0if**(X+1));
        xp = cexpf(*X + 1.0if**(X+1));
        xm = cexpf(-*X - 1.0if**(X+1));
        y = 0.5f * (xp-xm);
        //y = *X + 1.0if**(X+1); //this is probably ever so slightly faster but too close to tell
        //y = 0.5f * (cexpf(y)-cexpf(-y));
        *X = *(float *)&y; *++X = *((float *)&y+1);
    }
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int sinh_inplace_z (double *X, const size_t N)
{
    _Complex double y, xp, xm;
    
    for (size_t n=N; n>0u; --n, ++X)
    {
        xp = cexp(*X + 1.0i**(X+1));
        xm = cexp(-*X - 1.0i**(X+1));
        y = 0.5 * (xp-xm);
        *X = *(double *)&y; *++X = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
