//Gets natural logarithm (ln) of input X element-wise.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int log_s (float *Y, const float *X, const size_t N);
int log_d (double *Y, const double *X, const size_t N);
int log_c (float *Y, const float *X, const size_t N);
int log_z (double *Y, const double *X, const size_t N);

int log_inplace_s (float *X, const size_t N);
int log_inplace_d (double *X, const size_t N);
int log_inplace_c (float *X, const size_t N);
int log_inplace_z (double *X, const size_t N);


int log_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = logf(X[n]); }

    return 0;
}


int log_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { Y[n] = log(X[n]); }
    
    return 0;
}


int log_c (float *Y, const float *X, const size_t N)
{
    _Complex float y;

    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    for (size_t n2=0; n2<2*N; n2+=2)
    {
        y = clogf(X[n2]+1.0if*X[n2+1]);
        memcpy(&Y[n2],(float *)&y,2*sizeof(float));
        //Y[n2] = crealf(y); Y[n2+1] = cimagf(y);
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


int log_z (double *Y, const double *X, const size_t N)
{
    _Complex double y;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        y = clog(X[n2]+1.0i*X[n2+1]);
        memcpy(&Y[n2],(double *)&y,2*sizeof(double));
    }
    
    return 0;
}


int log_inplace_s (float *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = logf(X[n]); }

    return 0;
}


int log_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; n++) { X[n] = log(X[n]); }
    
    return 0;
}


int log_inplace_c (float *X, const size_t N)
{
    _Complex float x;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        x = clogf(X[n2]+1.0if*X[n2+1]);
        memcpy(&X[n2],(float *)&x,2*sizeof(float));
    }
    
    return 0;
}


int log_inplace_z (double *X, const size_t N)
{
    _Complex double x;

    for (size_t n2=0; n2<2*N; n2+=2)
    {
        x = clog(X[n2]+1.0i*X[n2+1]);
        memcpy(&X[n2],(double *)&x,2*sizeof(double));
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
