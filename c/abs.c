//This just gets absolute value of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <math.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int abs_s (float *Y, const float *X, const size_t N);
int abs_d (double *Y, const double *X, const size_t N);
int abs_c (float *Y, const float *X, const size_t N);
int abs_z (double *Y, const double *X, const size_t N);

int abs_inplace_s (float *X, const size_t N);
int abs_inplace_d (double *X, const size_t N);
int abs_inplace_c (float *X, const size_t N);
int abs_inplace_z (double *X, const size_t N);


int abs_s (float *Y, const float *X, const size_t N)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    for (size_t n=0; n<N; n++) { Y[n] = fabsf(X[n]); }
    //for (size_t n=0; n<N; n++) { Y[n] = (X[n]<0.0f) ? -X[n] : X[n]; }
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int abs_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = fabs(X[n]); }
    
    return 0;
}


int abs_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++, X+=2)
    {
        *Y++ = sqrtf(*X**X + *(X+1)**(X+1));
    }
    
    return 0;
}


int abs_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++, X+=2)
    {
        *Y++ = sqrt(*X**X + *(X+1)**(X+1));
    }
    
    return 0;
}


int abs_inplace_s (float *X, const size_t N)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    for (size_t n=0; n<N; n++) { X[n] = fabsf(X[n]); }
    //for (size_t n=0; n<N; n++) { if (X[n]<0.0f) { X[n] = -X[n]; } }
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int abs_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = fabs(X[n]); }
    
    return 0;
}


int abs_inplace_c (float *X, const size_t N)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    for (size_t n=0, n2=0; n<N; n++, n2+=2)
    {
        X[n] = sqrtf(X[n2]*X[n2] + X[n2+1]*X[n2+1]);
    }
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int abs_inplace_z (double *X, const size_t N)
{
    for (size_t n=0, n2=0; n<N; n++, n2+=2)
    {
        X[n] = sqrt(X[n2]*X[n2] + X[n2+1]*X[n2+1]);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
