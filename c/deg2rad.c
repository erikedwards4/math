//Converts degrees to radians for each element of X.
//This has in-place and not-in-place versions.

//LAPACKE_slascl was definitely slower than cblas_sscal.

#include <stdio.h>
#include <math.h>
//#include <lapacke.h>
//#include <time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int deg2rad_s (float *Y, const float *X, const size_t N);
int deg2rad_d (double *Y, const double *X, const size_t N);

int deg2rad_inplace_s (float *X, const size_t N);
int deg2rad_inplace_d (double *X, const size_t N);


int deg2rad_s (float *Y, const float *X, const size_t N)
{
    const float d2r = (float)(M_PI/180.0);

    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * d2r; }

    return 0;
}


int deg2rad_d (double *Y, const double *X, const size_t N)
{
    const double d2r = M_PI / 180.0;

    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * d2r; }
    
    return 0;
}



int deg2rad_inplace_s (float *X, const size_t N)
{
    const float d2r = (float)(M_PI/180.0);

    //if (LAPACKE_slascl_work(LAPACK_COL_MAJOR,'G',0,0,180.0f,(float)M_PI,1,N,X,1))
    //{ fprintf(stderr,"error in deg2rad_inplace_s: problem with LAPACKE function\n"); return 1; }

    //cblas_sscal((int)N,d2r,X,1);

    for (size_t n=N; n>0u; --n, ++X) { *X *= d2r; }

    return 0;
}


int deg2rad_inplace_d (double *X, const size_t N)
{
    const double d2r = M_PI / 180.0;

    for (size_t n=N; n>0u; --n, ++X) { *X *= d2r; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
