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
    if (N<1) { fprintf(stderr,"error in linspace_s: N (num elements Y) must be positive\n"); return 1; }
    
    const float stp = (b-a)/(N-1);

    Y[0] = a;
    for (int n=1; n<N-1; n++) { Y[n] = a + n*stp; }
    Y[N-1] = b;

    return 0;
}


int linspace_d (double *Y, const int N, const double a, const double b)
{
    if (N<1) { fprintf(stderr,"error in linspace_d: N (num elements Y) must be positive\n"); return 1; }

    const double stp = (b-a)/(N-1);

    Y[0] = a;
    for (int n=1; n<N-1; n++) { Y[n] = a + n*stp; }
    Y[N-1] = b;

    return 0;
}


int linspace_c (float *Y, const int N, const float a, const float b)
{
    if (N<1) { fprintf(stderr,"error in linspace_c: N (num elements Y) must be positive\n"); return 1; }

    const float z = 0.0f, stp = (b-a)/(N-1);

    Y[0] = a;
    for (int n=1; n<N-1; n++) { Y[2*n] = a + n*stp; }
    Y[2*(N-1)] = b;
    cblas_scopy(N,&z,0,&Y[1],2);  //set imag part to 0

    return 0;
}


int linspace_z (double *Y, const int N, const double a, const double b)
{
    if (N<1) { fprintf(stderr,"error in linspace_z: N (num elements Y) must be positive\n"); return 1; }

    const double z = 0.0, stp = (b-a)/(N-1);

    Y[0] = a;
    for (int n=1; n<N-1; n++) { Y[2*n] = a + n*stp; }
    Y[2*(N-1)] = b;
    cblas_dcopy(N,&z,0,&Y[1],2);  //set imag part to 0

    return 0;
}


#ifdef __cplusplus
}
}
#endif
