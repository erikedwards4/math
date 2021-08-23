//Gets conj part of complex-valued input X.
//This has in-place and not-in-place versions.

//LAPACKE_?lacgv was surprisingly slow!

#include <stdio.h>
//#include <lapacke.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int conj_c (float *Y, const float *X, const size_t N);
int conj_z (double *Y, const double *X, const size_t N);

int conj_inplace_c (float *X, const size_t N);
int conj_inplace_z (double *X, const size_t N);


int conj_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; *++Y = -*++X; }

    return 0;
}


int conj_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; *++Y = -*++X; }

    //cblas_dcopy(2*(int)N,X,1,Y,1);
    //cblas_dscal((int)N,-1.0,&Y[1],2);
    //for (n=0; n<2u*N; n+=2) { Y[n] = X[n]; }
    //for (n=1; n<2u*N; n+=2) { Y[n] = -X[n]; }

    return 0;
}


int conj_inplace_c (float *X, const size_t N)
{
    ++X;
    for (size_t n=N; n>0u; --n, X+=2) { *X = -*X; }
    
    return 0;
}


int conj_inplace_z (double *X, const size_t N)
{
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //LAPACKE_zlacgv_work((int)N,(_Complex double *)X,1);
    //cblas_dscal((int)N,-1.0,&X[1],2);
    ++X;
    for (size_t n=N; n>0u; --n, X+=2) { *X = -*X; }
    
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
