//Sets all N elements of Y equal to ipi.
//For complex cases, only real part is set to ipi.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ipi_s (float *Y, const int N);
int ipi_d (double *Y, const int N);
int ipi_c (float *Y, const int N);
int ipi_z (double *Y, const int N);


int ipi_s (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in ipi_s: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v = (float)M_1_PI;

    cblas_scopy(N,&v,0,Y,1);

    return 0;
}


int ipi_d (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in ipi_d: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v = M_1_PI;

    cblas_dcopy(N,&v,0,Y,1);

    return 0;
}


int ipi_c (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in ipi_c: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v[2] = {(float)M_1_PI,0.0f};

    cblas_ccopy(N,v,0,Y,1);

    return 0;
}


int ipi_z (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in ipi_z: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v[2] = {M_1_PI,0.0};

    cblas_zcopy(N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
