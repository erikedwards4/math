//Gets complex proj (projection to Riemann sphere) of complex-valued input X.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include "codee_math.h"
//#include <string.h>
#include <complex.h>
#include "codee_math.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int proj_c (float *Y, const float *X, const size_t N)
{
    _Complex float y;

    for (size_t n=N; n>0u; --n, X+=2, ++Y)
    {
        y = cprojf(*X + 1.0if**(X+1));
        *Y = *(float *)&y; *++Y = *((float *)&y+1);
    }

    return 0;
}


int proj_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=N; n>0u; --n, X+=2, ++Y)
    {
        y = cproj(*X + 1.0i**(X+1));
        *Y = *(double *)&y; *++Y = *((double *)&y+1);
    }
    
    return 0;
}


int proj_inplace_c (float *X, const size_t N)
{
    _Complex float y;

    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    //for (size_t n=0u; n<2u*N; n+=2)
    for (size_t n=N; n>0u; --n, ++X)
    {
        y = cprojf(*X + 1.0if**(X+1));
        //y = cprojf(X[n]+1.0if*X[n+1]);
        //memcpy(&X[n],(float *)&y,2*sizeof(float));
        //X[n] = *(float *)&y; X[n+1] = *((float *)&y+1);
        *X = *(float *)&y; *++X = *((float *)&y+1);
    }

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int proj_inplace_z (double *X, const size_t N)
{
    _Complex double y;

    for (size_t n=N; n>0u; --n, ++X)
    {
        y = cproj(*X + 1.0i**(X+1));
        *X = *(double *)&y; *++X = *((double *)&y+1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
