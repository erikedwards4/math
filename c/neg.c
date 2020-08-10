//Flips sign of input X element-wise: Y = -X.
//This has in-place and not-in-place versions.

//Using cblas_sscal was defintely slower than Y = -X,
//and in-place version was defintely faster.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int neg_s (float *Y, const float *X, const size_t N);
int neg_d (double *Y, const double *X, const size_t N);
int neg_c (float *Y, const float *X, const size_t N);
int neg_z (double *Y, const double *X, const size_t N);

int neg_inplace_s (float *X, const size_t N);
int neg_inplace_d (double *X, const size_t N);
int neg_inplace_c (float *X, const size_t N);
int neg_inplace_z (double *X, const size_t N);


int neg_s (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = -*X; }

    return 0;
}


int neg_d (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = -*X; }
    
    return 0;
}


int neg_c (float *Y, const float *X, const size_t N)
{
    for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = -*X; }
    
    return 0;
}


int neg_z (double *Y, const double *X, const size_t N)
{
    for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = -*X; }
    
    return 0;
}


int neg_inplace_s (float *X, const size_t N)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    for (size_t n=0; n<N; ++n, ++X) { *X = -*X; }
    //cblas_sscal((int)N,-1.0f,X,1);
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int neg_inplace_d (double *X, const size_t N)
{
    for (size_t n=0; n<N; ++n, ++X) { *X = -*X; }
    
    return 0;
}


int neg_inplace_c (float *X, const size_t N)
{
    for (size_t n=0; n<2*N; ++n, ++X) { *X = -*X; }
    
    return 0;
}


int neg_inplace_z (double *X, const size_t N)
{
    for (size_t n=0; n<2*N; ++n, ++X) { *X = -*X; }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
