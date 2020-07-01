//Sets all N elements of Y equal to e.
//For complex cases, only real part is set to e.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int e_s (float *Y, const int N);
int e_d (double *Y, const int N);
int e_c (float *Y, const int N);
int e_z (double *Y, const int N);


int e_s (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in e_s: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v = (float)M_E;

    cblas_scopy(N,&v,0,Y,1);

    return 0;
}


int e_d (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in e_d: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v = M_E;

    cblas_dcopy(N,&v,0,Y,1);

    return 0;
}


int e_c (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in e_c: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v[2] = {(float)M_E,0.0f};

    cblas_ccopy(N,v,0,Y,1);

    return 0;
}


int e_z (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in e_z: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v[2] = {M_E,0.0};

    cblas_zcopy(N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
