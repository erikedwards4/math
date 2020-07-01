//Generates vector Y with elements linearly spaced from a to b.
//For complex Y, imag part is set to 0.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int linspace_s (float *Y, const int N, const float a, const float b);
int linspace_d (double *Y, const int N, const double a, const double b);
int linspace_c (float *Y, const int N, const float a, const float b);
int linspace_z (double *Y, const int N, const double a, const double b);


int linspace_s (float *Y, const int N, const float a, const float b)
{
    const float stp = (b-a)/(N-1);
    int n;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<1) { fprintf(stderr,"error in linspace_s: N (num elements Y) must be positive\n"); return 1; }

    Y[0] = a;
    for (n=1; n<N-1; n++) { Y[n] = a + n*stp; }
    Y[N-1] = b;

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int linspace_d (double *Y, const int N, const double a, const double b)
{
    const double stp = (b-a)/(N-1);
    int n;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<1) { fprintf(stderr,"error in linspace_d: N (num elements Y) must be positive\n"); return 1; }

    Y[0] = a;
    for (n=1; n<N-1; n++) { Y[n] = a + n*stp; }
    Y[N-1] = b;

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int linspace_c (float *Y, const int N, const float a, const float b)
{
    const float z = 0.0f, stp = (b-a)/(N-1);
    int n;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<1) { fprintf(stderr,"error in linspace_c: N (num elements Y) must be positive\n"); return 1; }

    //Set imag part to 0
    cblas_scopy(N,&z,0,&Y[1],2);

    //Set real part
    Y[0] = a;
    for (n=1; n<N-1; n++) { Y[2*n] = a + n*stp; }
    Y[2*(N-1)] = b;

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int linspace_z (double *Y, const int N, const double a, const double b)
{
    const double z = 0.0, stp = (b-a)/(N-1);
    int n;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (N<1) { fprintf(stderr,"error in linspace_z: N (num elements Y) must be positive\n"); return 1; }

    //Set imag part to 0
    cblas_dcopy(N,&z,0,&Y[1],2);

    //Set real part
    Y[0] = a;
    for (n=1; n<N-1; n++) { Y[2*n] = a + n*stp; }
    Y[2*(N-1)] = b;

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
