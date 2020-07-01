//Sets all elements of Y equal to 0, except 1s on the main diagonal.
//LAPACKE_slaset was definitely slower than the cblas_scopy solution. 

#include <stdio.h>
#include <cblas.h>
//#include <lapacke.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int eye_s (float *Y, const int R, const int C, const char iscolmajor);
int eye_d (double *Y, const int R, const int C, const char iscolmajor);
int eye_c (float *Y, const int R, const int C, const char iscolmajor);
int eye_z (double *Y, const int R, const int C, const char iscolmajor);


int eye_s (float *Y, const int R, const int C, const char iscolmajor)
{
    if (R<0) { fprintf(stderr,"error in eye_s: R (num rows Y) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in eye_s: C (num cols Y) must be nonnegative\n"); return 1; }

    const float z = 0.0f, o = 1.0f;
    const int M = (R<C) ? R : C;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    cblas_scopy(R*C,&z,0,Y,1);
    if (iscolmajor) { cblas_scopy(M,&o,0,Y,R+1); }
    else { cblas_scopy(M,&o,0,Y,C+1); }

    //LAPACKE solution
    //const int LO = (iscolmajor) ? LAPACK_COL_MAJOR : LAPACK_ROW_MAJOR;
    //const int lda = (iscolmajor) ? R : C;
    //if (LAPACKE_slaset_work(LO,'A',R,C,0.0f,1.0f,Y,lda))
    //{ fprintf(stderr,"error in eye_s: problem with LAPACKE function\n"); return 1; }

    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int eye_d (double *Y, const int R, const int C, const char iscolmajor)
{
    if (R<0) { fprintf(stderr,"error in eye_d: R (num rows Y) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in eye_d: C (num cols Y) must be nonnegative\n"); return 1; }

    const double z = 0.0, o = 1.0;
    const int M = (R<C) ? R : C;

    cblas_dcopy(R*C,&z,0,Y,1);
    if (iscolmajor) { cblas_dcopy(M,&o,0,Y,R+1); }
    else { cblas_dcopy(M,&o,0,Y,C+1); }

    return 0;
}


int eye_c (float *Y, const int R, const int C, const char iscolmajor)
{
    if (R<0) { fprintf(stderr,"error in eye_c: R (num rows Y) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in eye_c: C (num cols Y) must be nonnegative\n"); return 1; }

    const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
    const int M = (R<C) ? R : C;

    cblas_ccopy(R*C,z,0,Y,1);
    if (iscolmajor) { cblas_ccopy(M,o,0,Y,R+1); }
    else { cblas_ccopy(M,o,0,Y,C+1); }

    return 0;
}


int eye_z (double *Y, const int R, const int C, const char iscolmajor)
{
    if (R<0) { fprintf(stderr,"error in eye_z: R (num rows Y) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in eye_z: C (num cols Y) must be nonnegative\n"); return 1; }

    const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
    const int M = (R<C) ? R : C;

    cblas_zcopy(R*C,z,0,Y,1);
    if (iscolmajor) { cblas_zcopy(M,o,0,Y,R+1); }
    else { cblas_zcopy(M,o,0,Y,C+1); }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
