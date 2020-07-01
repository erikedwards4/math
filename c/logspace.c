//Generates vector Y with elements logarithmically spaced from 10^a to 10^b.
//For complex Y, imag part is set to 0.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int logspace_s (float *Y, const int N, const float a, const float b);
int logspace_d (double *Y, const int N, const double a, const double b);
int logspace_c (float *Y, const int N, const float a, const float b);
int logspace_z (double *Y, const int N, const double a, const double b);


int logspace_s (float *Y, const int N, const float a, const float b)
{
    if (N<1) { fprintf(stderr,"error in logspace_s: N (num elements Y) must be positive\n"); return 1; }
    
    const float stp = (b-a)/(N-1);
    int n;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    Y[0] = a;
    for (n=1; n<N-1; n++) { Y[n] = a + n*stp; }
    Y[N-1] = b;
    for (n=0; n<N; n++) { Y[n] = powf(10.0f,Y[n]); }

    // Y[0] = powf(10.0f,a);
    // for (n=1; n<N-1; n++) { Y[n] = powf(10.0f,a+n*stp); }
    // Y[N-1] = powf(10.0f,b);

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int logspace_d (double *Y, const int N, const double a, const double b)
{
    if (N<1) { fprintf(stderr,"error in logspace_d: N (num elements Y) must be positive\n"); return 1; }

    const double stp = (b-a)/(N-1);
    int n;

    Y[0] = a;
    for (n=1; n<N-1; n++) { Y[n] = a + n*stp; }
    Y[N-1] = b;
    for (n=0; n<N; n++) { Y[n] = pow(10.0,Y[n]); }

    return 0;
}


int logspace_c (float *Y, const int N, const float a, const float b)
{
    if (N<1) { fprintf(stderr,"error in logspace_c: N (num elements Y) must be positive\n"); return 1; }
    const float z = 0.0f, stp = (b-a)/(N-1);
    int n;

    Y[0] = a;
    for (n=1; n<N-1; n++) { Y[2*n] = a + n*stp; }
    Y[N-1] = b;
    for (n=0; n<N; n++) { Y[2*n] = powf(10.0f,Y[2*n]); }
    cblas_scopy(N,&z,0,&Y[1],2);  //set imag part to 0

    return 0;
}


int logspace_z (double *Y, const int N, const double a, const double b)
{
    if (N<1) { fprintf(stderr,"error in logspace_z: N (num elements Y) must be positive\n"); return 1; }

    const double z = 0.0, stp = (b-a)/(N-1);
    int n;

    Y[0] = a;
    for (n=1; n<N-1; n++) { Y[2*n] = a + n*stp; }
    Y[N-1] = b;
    for (n=0; n<N; n++) { Y[2*n] = pow(10.0,Y[2*n]); }
    cblas_dcopy(N,&z,0,&Y[1],2);  //set imag part to 0

    return 0;
}


#ifdef __cplusplus
}
}
#endif
