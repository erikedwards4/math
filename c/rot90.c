//Rotates matrix X by 90 degrees K times.


#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rot90_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int K);
int rot90_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int K);
int rot90_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int K);
int rot90_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int K);


int rot90_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int K)
{
    const size_t N = R*C;

    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    if (N==0) {}
    else if (K%4==0 || R==1 || C==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (K%4==1)
    {
        if (iscolmajor)
        {}
        else
        {}
    }
    else if (K%4==2)
    {
        if (iscolmajor)
        {
            Y += N - 1;
            for (n=0; n<N; ++n, ++X, --Y) { *Y = *X; }
        }
        else
        {}
    }
    else // (K%4==3)
    {
        if (iscolmajor)
        {}
        else
        {}
    }
    
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int rot90_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int K)
{
    if (R==1 || C==1) { cblas_dcopy((int)(R*C),X,1,Y,1); }
    else if (iscolmajor)
    {
        for (size_t c=0; c<C; ++c, X+=R, ++Y) { cblas_dcopy((int)(R),X,1,Y,(int)C); }
    }
    else
    {
        for (size_t r=0; r<R; ++r, X+=C, ++Y) { cblas_dcopy((int)(C),X,1,Y,(int)R); }
    }

    return 0;
}


int rot90_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor, const int K)
{
    if (R==1 || C==1) { cblas_ccopy((int)(R*C),X,1,Y,1); }
    else if (iscolmajor)
    {
        for (size_t c=0; c<C; c+=2, X+=2*R, Y+=2) { cblas_ccopy((int)(R),X,1,Y,(int)C); }
    }
    else
    {
        for (size_t r=0; r<R; r+=2, X+=2*C, Y+=2) { cblas_ccopy((int)(C),X,1,Y,(int)R); }
    }

    return 0;
}


int rot90_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor, const int K)
{
    if (R==1 || C==1) { cblas_zcopy((int)(R*C),X,1,Y,1); }
    else if (iscolmajor)
    {
        for (size_t c=0; c<C; c+=2, X+=2*R, Y+=2) { cblas_zcopy((int)(R),X,1,Y,(int)C); }
    }
    else
    {
        for (size_t r=0; r<R; r+=2, X+=2*C, Y+=2) { cblas_zcopy((int)(C),X,1,Y,(int)R); }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
