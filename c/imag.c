//Gets imaginary part of complex-valued input X (output is real-valued).

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int imag_c (float *Y, const float *X, const size_t N);
int imag_z (double *Y, const double *X, const size_t N);


int imag_c (float *Y, const float *X, const size_t N)
{
    cblas_scopy((int)N,&X[1],2,Y,1);

    return 0;
}


int imag_z (double *Y, const double *X, const size_t N)
{
    cblas_dcopy((int)N,&X[1],2,Y,1);
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
