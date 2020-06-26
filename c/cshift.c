//Circular shift of X by P elements along dim.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cshift_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int cshift_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int cshift_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int cshift_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);

int cshift_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int cshift_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int cshift_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int cshift_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);


int cshift_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    int l, n = 0;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (R<0) { fprintf(stderr,"error in cshift_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in cshift_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in cshift_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in cshift_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * I;
    const int L = R*C*S*H/(I*N);
    const int D = P % N;

    //Sort
    if (D==0) { cblas_scopy(R*C*S*H,X,1,Y,1); }
    else if (D<0)
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_scopy((N+D)*I,&X[n-D*I],1,&Y[n],1);
            cblas_scopy(-D*I,&X[n],1,&Y[n+(N+D)*I],1);
        }
    }
    else
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_scopy((N-D)*I,&X[n],1,&Y[n+D*I],1);
            cblas_scopy(D*I,&X[n+(N-D)*I],1,&Y[n],1);
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int cshift_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    int l, n = 0;

    //Checks
    if (R<0) { fprintf(stderr,"error in cshift_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in cshift_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in cshift_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in cshift_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * I;
    const int L = R*C*S*H/(I*N);
    const int D = P % N;

    //Sort
    if (D==0) { cblas_dcopy(R*C*S*H,X,1,Y,1); }
    else if (D<0)
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_dcopy((N+D)*I,&X[n-D*I],1,&Y[n],1);
            cblas_dcopy(-D*I,&X[n],1,&Y[n+(N+D)*I],1);
        }
    }
    else
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_dcopy((N-D)*I,&X[n],1,&Y[n+D*I],1);
            cblas_dcopy(D*I,&X[n+(N-D)*I],1,&Y[n],1);
        }
    }

    return 0;
}


int cshift_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    int l, n = 0;

    //Checks
    if (R<0) { fprintf(stderr,"error in cshift_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in cshift_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in cshift_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in cshift_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * I;
    const int L = R*C*S*H/(I*N);
    const int D = P % N;

    //Sort
    if (D==0) { cblas_ccopy(R*C*S*H,X,1,Y,1); }
    else if (D<0)
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((N+D)*I,&X[n-2*D*I],1,&Y[n],1);
            cblas_ccopy(-D*I,&X[n],1,&Y[n+2*(N+D)*I],1);
        }
    }
    else
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((N-D)*I,&X[n],1,&Y[n+2*D*I],1);
            cblas_ccopy(D*I,&X[n+2*(N-D)*I],1,&Y[n],1);
        }
    }

    return 0;
}


int cshift_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    int l, n = 0;

    //Checks
    if (R<0) { fprintf(stderr,"error in cshift_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in cshift_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in cshift_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in cshift_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * I;
    const int L = R*C*S*H/(I*N);
    const int D = P % N;

    //Sort
    if (D==0) { cblas_zcopy(R*C*S*H,X,1,Y,1); }
    else if (D<0)
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((N+D)*I,&X[n-2*D*I],1,&Y[n],1);
            cblas_zcopy(-D*I,&X[n],1,&Y[n+2*(N+D)*I],1);
        }
    }
    else
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((N-D)*I,&X[n],1,&Y[n+2*D*I],1);
            cblas_zcopy(D*I,&X[n+2*(N-D)*I],1,&Y[n],1);
        }
    }

    return 0;
}


int cshift_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    int l, n = 0;
    float *X1;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (R<0) { fprintf(stderr,"error in cshift_inplace_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in cshift_inplace_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in cshift_inplace_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in cshift_inplace_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * I;
    const int L = R*C*S*H/(I*N);
    const int D = P % N;

    //Allocate
    const size_t MS = (D>0) ? (size_t)(I*(N-D)) : (size_t)(-D*I);
    if (!(X1=(float *)malloc(MS*sizeof(float)))) { fprintf(stderr,"error in cshift_inplace_s: problem with malloc. "); perror("malloc"); return 1; }

    //Sort
    if (D<0)
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_scopy(-D*I,&X[n],1,X1,1);
            cblas_scopy((N+D)*I,&X[n-D*I],1,&X[n],1);
            cblas_scopy(-D*I,X1,1,&X[n+(N+D)*I],1);
        }
    }
    else if (D>0)
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_scopy((N-D)*I,&X[n],1,X1,1);
            cblas_scopy(D*I,&X[n+(N-D)*I],1,&X[n],1);
            cblas_scopy((N-D)*I,X1,1,&X[n+D*I],1);
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    free(X1);
    return 0;
}


int cshift_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    int l, n = 0;
    double *X1;

    //Checks
    if (R<0) { fprintf(stderr,"error in cshift_inplace_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in cshift_inplace_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in cshift_inplace_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in cshift_inplace_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * I;
    const int L = R*C*S*H/(I*N);
    const int D = P % N;

    //Allocate
    const size_t MS = (D>0) ? (size_t)(I*(N-D)) : (size_t)(-D*I);
    if (!(X1=(double *)malloc(MS*sizeof(double)))) { fprintf(stderr,"error in cshift_inplace_d: problem with malloc. "); perror("malloc"); return 1; }

    //Sort
    if (D<0)
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_dcopy(-D*I,&X[n],1,X1,1);
            cblas_dcopy((N+D)*I,&X[n-D*I],1,&X[n],1);
            cblas_dcopy(-D*I,X1,1,&X[n+(N+D)*I],1);
        }
    }
    else if (D>0)
    {
        for (l=0; l<L; l++, n+=J)
        {
            cblas_dcopy((N-D)*I,&X[n],1,X1,1);
            cblas_dcopy(D*I,&X[n+(N-D)*I],1,&X[n],1);
            cblas_dcopy((N-D)*I,X1,1,&X[n+D*I],1);
        }
    }

    free(X1);
    return 0;
}


int cshift_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    int l, n = 0;
    float *X1;

    //Checks
    if (R<0) { fprintf(stderr,"error in cshift_inplace_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in cshift_inplace_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in cshift_inplace_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in cshift_inplace_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * I;
    const int L = R*C*S*H/(I*N);
    const int D = P % N;

    //Allocate
    const size_t MS = (D>0) ? (size_t)(I*(N-D)) : (size_t)(-D*I);
    if (!(X1=(float *)malloc(MS*2u*sizeof(float)))) { fprintf(stderr,"error in cshift_inplace_c: problem with malloc. "); perror("malloc"); return 1; }

    //Sort
    if (D<0)
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy(-D*I,&X[n],1,X1,1);
            cblas_ccopy((N+D)*I,&X[n-2*D*I],1,&X[n],1);
            cblas_ccopy(-D*I,X1,1,&X[n+2*(N+D)*I],1);
        }
    }
    else if (D>0)
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((N-D)*I,&X[n],1,X1,1);
            cblas_ccopy(D*I,&X[n+2*(N-D)*I],1,&X[n],1);
            cblas_ccopy((N-D)*I,X1,1,&X[n+2*D*I],1);
        }
    }

    free(X1);
    return 0;
}


int cshift_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    int l, n = 0;
    double *X1;

    //Checks
    if (R<0) { fprintf(stderr,"error in cshift_inplace_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in cshift_inplace_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in cshift_inplace_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in cshift_inplace_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * I;
    const int L = R*C*S*H/(I*N);
    const int D = P % N;

    //Allocate
    const size_t MS = (D>0) ? (size_t)(I*(N-D)) : (size_t)(-D*I);
    if (!(X1=(double *)malloc(MS*2u*sizeof(double)))) { fprintf(stderr,"error in cshift_inplace_z: problem with malloc. "); perror("malloc"); return 1; }

    //Sort
    if (D<0)
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy(-D*I,&X[n],1,X1,1);
            cblas_zcopy((N+D)*I,&X[n-2*D*I],1,&X[n],1);
            cblas_zcopy(-D*I,X1,1,&X[n+2*(N+D)*I],1);
        }
    }
    else if (D>0)
    {
        for (l=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((N-D)*I,&X[n],1,X1,1);
            cblas_zcopy(D*I,&X[n+2*(N-D)*I],1,&X[n],1);
            cblas_zcopy((N-D)*I,X1,1,&X[n+2*D*I],1);
        }
    }

    free(X1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
