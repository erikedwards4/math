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
int square_inplace_c (float *X, const size_t N);
int square_inplace_z (double *X, const size_t N);


int square_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = X[n]*X[n]; }

    return 0;
}


int square_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = X[n]*X[n]; }
    
    return 0;
}


int square_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++, X+=2) { *Y++ = *X**X + *(X+1)**(X+1); }
    
    return 0;
}


int square_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++, X+=2) { *Y++ = *X**X + *(X+1)**(X+1); }
    
    return 0;
}


int square_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] *= X[n]; }

    return 0;
}


int square_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] *= X[n]; }
    
    return 0;
}


int square_inplace_c (float *X, const size_t N)
{
    for (size_t n=0, n2=0; n<N; n++, n2+=2) { X[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1]; }
    
    return 0;
}


int square_inplace_z (double *X, const size_t N)
{
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    //for (size_t n=0; n<N; n++) { X[n] = X[n*2]*X[n*2] + X[n*2+1]*X[n*2+1]; }
    for (size_t n=0, n2=0; n<N; n++, n2+=2) { X[n] = X[n2]*X[n2] + X[n2+1]*X[n2+1]; }
    //for (size_t n=0; n<2*N; n++) { X[n] *= X[n]; }
    //for (size_t n=0; n<N; n++) { X[n] = X[2*n] + X[2*n+1]; }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
