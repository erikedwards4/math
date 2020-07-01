//Sets all N elements of Y equal to realmax (FLT_MAX or DBL_MAX).
//For complex cases, only real part is set to realmax.

#include <stdio.h>
#include <float.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int realmax_s (float *Y, const int N);
int realmax_d (double *Y, const int N);
int realmax_c (float *Y, const int N);
int realmax_z (double *Y, const int N);


int realmax_s (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in realmax_s: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v = FLT_MAX;

    cblas_scopy(N,&v,0,Y,1);

    return 0;
}


int realmax_d (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in realmax_d: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v = DBL_MAX;

    cblas_dcopy(N,&v,0,Y,1);

    return 0;
}


int realmax_c (float *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in realmax_c: N (num elements Y) must be nonnegative\n"); return 1; }

    const float v[2] = {FLT_MAX,0.0f};

    cblas_ccopy(N,v,0,Y,1);

    return 0;
}


int realmax_z (double *Y, const int N)
{
    if (N<0) { fprintf(stderr,"error in realmax_z: N (num elements Y) must be nonnegative\n"); return 1; }

    const double v[2] = {DBL_MAX,0.0};

    cblas_zcopy(N,v,0,Y,1);

    return 0;
}


#ifdef __cplusplus
}
}
#endif
