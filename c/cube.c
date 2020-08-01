//This just cubes input X element-wise.
//For complex data, this is |X|.^2, i.e. Xr*Xr + Xi*Xi.

//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cube_s (float *Y, const float *X, const size_t N);
int cube_d (double *Y, const double *X, const size_t N);
int cube_c (float *Y, const float *X, const size_t N);
int cube_z (double *Y, const double *X, const size_t N);

int cube_inplace_s (float *X, const size_t N);
int cube_inplace_d (double *X, const size_t N);
int cube_inplace_c (float *X, const size_t N);
int cube_inplace_z (double *X, const size_t N);


int cube_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X * *X * *X; }

    return 0;
}


int cube_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X * *X * *X; }
    
    return 0;
}


int cube_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0, n2=0; n<N; ++n, n2+=2)
    {
        *Y = X[n2]*X[n2] + X[n2+1]*X[n2+1];
        *Y *= sqrtf(*Y);
    }
    
    return 0;
}


int cube_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0, n2=0; n<N; ++n, n2+=2)
    {
        *Y = X[n2]*X[n2] + X[n2+1]*X[n2+1];
        *Y *= sqrt(*Y);
    }
    
    return 0;
}


int cube_inplace_s (float *X, const size_t N)
{
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    for (size_t n=0; n<N; ++n, ++X) { *X = *X * *X * *X; }
    //for (size_t n=0; n<N; ++n, ++X) { *X *= *X * *X; }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int cube_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = *X * *X * *X; }
    
    return 0;
}


int cube_inplace_c (float *X, const size_t N)
{
    for (size_t n=0, n2=0; n<N; ++n, n2+=2)
    {
        *X = X[n2]*X[n2] + X[n2+1]*X[n2+1];
        *X *= sqrtf(*X);
    }
    
    return 0;
}


int cube_inplace_z (double *X, const size_t N)
{
    for (size_t n=0, n2=0; n<N; ++n, n2+=2)
    {
        *X = X[n2]*X[n2] + X[n2+1]*X[n2+1];
        *X *= sqrt(*X);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
