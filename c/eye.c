//Sets all elements of Y equal to 0, except 1s on the main diagonal.
//LAPACKE_slaset was definitely slower than the cblas_scopy solution.
//But, as usual, the direct for loop is definitely faster for small N,
//and slightly faster or same speed for large N.

#include <stdio.h>
//#include <cblas.h>
//#include <lapacke.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int eye_s (float *Y, const size_t R, const size_t C, const char iscolmajor);
int eye_d (double *Y, const size_t R, const size_t C, const char iscolmajor);
int eye_c (float *Y, const size_t R, const size_t C, const char iscolmajor);
int eye_z (double *Y, const size_t R, const size_t C, const char iscolmajor);


int eye_s (float *Y, const size_t R, const size_t C, const char iscolmajor)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    //const float z = 0.0f, o = 1.0f;
    //const size_t N = R*C, M = (R<C) ? R : C;

    const size_t D0 = (iscolmajor) ? C : R, D1 = (iscolmajor) ? R : C;
    for (size_t d0=0; d0<D0; ++d0)
    {
        for (size_t d1=0; d1<D1; ++d1, ++Y) { *Y = (d0==d1); }
    }

    // if (iscolmajor)
    // {
    //     for (size_t c=0; c<C; ++c)
    //     {
    //         for (size_t r=0; r<R; ++r, ++Y) { *Y = (r==c); }
    //     }
    // }
    // else
    // {
    //     for (size_t r=0; r<R; ++r)
    //     {
    //         for (size_t c=0; c<C; ++c, ++Y) { *Y = (r==c); }
    //     }
    // }

    // cblas_scopy((int)N,&z,0,Y,1);
    // if (iscolmajor) { cblas_scopy((int)M,&o,0,Y,(int)R+1); }
    // else { cblas_scopy((int)M,&o,0,Y,(int)C+1); }

    //LAPACKE solution
    //const size_t LO = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
    //const int lda = (iscolmajor) ? R : C;
    //if (LAPACKE_slaset_work(LO,'A',(int)R,(int)C,0.0f,1.0f,Y,lda))
    //{ fprintf(stderr,"error in eye_s: problem with LAPACKE function\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int eye_d (double *Y, const size_t R, const size_t C, const char iscolmajor)
{
    const size_t D0 = (iscolmajor) ? C : R, D1 = (iscolmajor) ? R : C;

    for (size_t d0=0; d0<D0; ++d0)
    {
        for (size_t d1=0; d1<D1; ++d1, ++Y) { *Y = (d0==d1); }
    }

    return 0;
}


int eye_c (float *Y, const size_t R, const size_t C, const char iscolmajor)
{
    const size_t D0 = (iscolmajor) ? C : R, D1 = (iscolmajor) ? R : C;
    
    for (size_t d0=0; d0<D0; ++d0)
    {
        for (size_t d1=0; d1<D1; ++d1, ++Y) { *Y = (d0==d1); *++Y = 0.0f; }
    }

    return 0;
}


int eye_z (double *Y, const size_t R, const size_t C, const char iscolmajor)
{
    const size_t D0 = (iscolmajor) ? C : R, D1 = (iscolmajor) ? R : C;
    
    for (size_t d0=0; d0<D0; ++d0)
    {
        for (size_t d1=0; d1<D1; ++d1, ++Y) { *Y = (d0==d1); *++Y = 0.0; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
