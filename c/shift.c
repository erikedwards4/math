//Shifts X by P elements along dim.
//The P edge samples are replaced by 0s.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int shift_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int shift_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int shift_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int shift_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);

int shift_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int shift_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int shift_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int shift_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);


int shift_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    const float z = 0.0f;
    int l, n = 0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (R<0) { fprintf(stderr,"error in shift_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in shift_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in shift_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in shift_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * K;
    const int L = R*C*S*H/(K*N);

    //Sort
    if (P==0) { cblas_scopy(R*C*S*H,X,1,Y,1); }
    else if (P<=-N || P>=N) { cblas_scopy(R*C*S*H,&z,0,Y,1); }
    else if (P<0)
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_scopy((N+P)*K,&X[n-P*K],1,&Y[n],1);
            cblas_scopy(-P*K,&z,0,&Y[n+(N+P)*K],1);
        }
    }
    else
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_scopy((N-P)*K,&X[n],1,&Y[n+P*K],1);
            cblas_scopy(P*K,&z,0,&Y[n],1);
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int shift_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    const double z = 0.0;
    int l, n = 0;

    //Checks
    if (R<0) { fprintf(stderr,"error in shift_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in shift_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in shift_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in shift_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * K;
    const int L = R*C*S*H/(K*N);

    //Sort
    if (P==0) { cblas_dcopy(R*C*S*H,X,1,Y,1); }
    else if (P<=-N || P>=N) { cblas_dcopy(R*C*S*H,&z,0,Y,1); }
    else if (P<0)
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_dcopy((N+P)*K,&X[n-P*K],1,&Y[n],1);
            cblas_dcopy(-P*K,&z,0,&Y[n+(N+P)*K],1);
        }
    }
    else
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_dcopy((N-P)*K,&X[n],1,&Y[n+P*K],1);
            cblas_dcopy(P*K,&z,0,&Y[n],1);
        }
    }
    
    return 0;
}


int shift_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    const float z[2] = {0.0f,0.0f};
    int l, n = 0;

    //Checks
    if (R<0) { fprintf(stderr,"error in shift_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in shift_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in shift_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in shift_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * K;
    const int L = R*C*S*H/(K*N);

    //Sort
    if (P==0) { cblas_ccopy(R*C*S*H,X,1,Y,1); }
    else if (P<=-N || P>=N) { cblas_ccopy(R*C*S*H,z,0,Y,1); }
    else if (P<0)
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((N+P)*K,&X[n-2*P*K],1,&Y[n],1);
            cblas_ccopy(-P*K,z,0,&Y[n+2*(N+P)*K],1);
        }
    }
    else
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((N-P)*K,&X[n],1,&Y[n+2*P*K],1);
            cblas_ccopy(P*K,z,0,&Y[n],1);
        }
    }
    
    return 0;
}


int shift_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    const double z[2] = {0.0,0.0};
    int l, n = 0;

    //Checks
    if (R<0) { fprintf(stderr,"error in shift_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in shift_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in shift_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in shift_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * K;
    const int L = R*C*S*H/(K*N);

    //Sort
    if (P==0) { cblas_zcopy(R*C*S*H,X,1,Y,1); }
    else if (P<=-N || P>=N) { cblas_zcopy(R*C*S*H,z,0,Y,1); }
    else if (P<0)
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((N+P)*K,&X[n-2*P*K],1,&Y[n],1);
            cblas_zcopy(-P*K,z,0,&Y[n+2*(N+P)*K],1);
        }
    }
    else
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((N-P)*K,&X[n],1,&Y[n+2*P*K],1);
            cblas_zcopy(P*K,z,0,&Y[n],1);
        }
    }
    
    return 0;
}


int shift_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    const float z = 0.0f;
    int l, n = 0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (R<0) { fprintf(stderr,"error in shift_inplace_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in shift_inplace_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in shift_inplace_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in shift_inplace_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * K;
    const int L = R*C*S*H/(K*N);

    //Sort
    if (P<=-N || P>=N) { cblas_scopy(R*C*S*H,&z,0,X,1); }
    else if (P<0)
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_scopy((N+P)*K,&X[n-P*K],1,&X[n],1);
            cblas_scopy(-P*K,&z,0,&X[n+(N+P)*K],1);
        }
    }
    else if (P>0)
    {
        float *X1;
        if (!(X1=(float *)malloc((size_t)(K*(N-P))*sizeof(float)))) { fprintf(stderr,"error in shift_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        for (l=0; l<L; l++, n+=J)
        {
            cblas_scopy((N-P)*K,&X[n],1,X1,1);
            cblas_scopy((N-P)*K,X1,1,&X[n+P*K],1);
            cblas_scopy(P*K,&z,0,&X[n],1);
        }
        free(X1);
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int shift_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    const double z = 0.0;
    int l, n = 0;

    //Checks
    if (R<0) { fprintf(stderr,"error in shift_inplace_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in shift_inplace_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in shift_inplace_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in shift_inplace_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * K;
    const int L = R*C*S*H/(K*N);

    //Sort
    if (P<=-N || P>=N) { cblas_dcopy(R*C*S*H,&z,0,X,1); }
    else if (P<0)
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_dcopy((N+P)*K,&X[n-P*K],1,&X[n],1);
            cblas_dcopy(-P*K,&z,0,&X[n+(N+P)*K],1);
        }
    }
    else if (P>0)
    {
        double *X1;
        if (!(X1=(double *)malloc((size_t)(K*(N-P))*sizeof(double)))) { fprintf(stderr,"error in shift_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        for (l=0; l<L; l++, n+=J)
        {
            cblas_dcopy((N-P)*K,&X[n],1,X1,1);
            cblas_dcopy((N-P)*K,X1,1,&X[n+P*K],1);
            cblas_dcopy(P*K,&z,0,&X[n],1);
        }
        free(X1);
    }

    return 0;
}


int shift_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    const float z[2] = {0.0f,0.0f};
    int l, n = 0;

    //Checks
    if (R<0) { fprintf(stderr,"error in shift_inplace_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in shift_inplace_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in shift_inplace_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in shift_inplace_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * K;
    const int L = R*C*S*H/(K*N);

    //Sort
    if (P<=-N || P>=N) { cblas_ccopy(R*C*S*H,z,0,X,1); }
    else if (P<0)
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((N+P)*K,&X[n-2*P*K],1,&X[n],1);
            cblas_ccopy(-P*K,z,0,&X[n+2*(N+P)*K],1);
        }
    }
    else if (P>0)
    {
        float *X1;
        if (!(X1=(float *)malloc((size_t)(2*K*(N-P))*sizeof(float)))) { fprintf(stderr,"error in shift_inplace_c: problem with malloc. "); perror("malloc"); return 1; }
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((N-P)*K,&X[n],1,X1,1);
            cblas_ccopy((N-P)*K,X1,1,&X[n+2*P*K],1);
            cblas_ccopy(P*K,z,0,&X[n],1);
        }
        free(X1);
    }
    
    return 0;
}


int shift_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    const double z[2] = {0.0,0.0};
    int l, n = 0;

    //Checks
    if (R<0) { fprintf(stderr,"error in shift_inplace_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in shift_inplace_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in shift_inplace_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in shift_inplace_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * K;
    const int L = R*C*S*H/(K*N);

    //Sort
    if (P<=-N || P>=N) { cblas_zcopy(R*C*S*H,z,0,X,1); }
    else if (P<0)
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((N+P)*K,&X[n-2*P*K],1,&X[n],1);
            cblas_zcopy(-P*K,z,0,&X[n+2*(N+P)*K],1);
        }
    }
    else if (P>0)
    {
        double *X1;
        if (!(X1=(double *)malloc((size_t)(2*K*(N-P))*sizeof(double)))) { fprintf(stderr,"error in shift_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((N-P)*K,&X[n],1,X1,1);
            cblas_zcopy((N-P)*K,X1,1,&X[n+2*P*K],1);
            cblas_zcopy(P*K,z,0,&X[n],1);
        }
        free(X1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
