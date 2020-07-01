//Prepads X with P zeros along dim.
//Thus, Y has the same size as X along all other dims,
//but has a greater length than X along dim.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <cblas.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int prepad_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int prepad_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int prepad_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);
int prepad_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P);


int prepad_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    if (R<0) { fprintf(stderr,"error in prepad_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in prepad_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in prepad_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in prepad_s: H (num hyperslices X) must be nonnegative\n"); return 1; }
    if (P<0) { fprintf(stderr,"error in prepad_s: P (num prepads) must be nonnegative\n"); return 1; }

    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int N2 = N1 + P;
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J1 = N1*I, J2 = N2*I;
    const int L1 = R*C*S*H/(I*N1);
    const int L2 = (dim==0) ? (R+P)*C*S*H/(I*N2) : (dim==1) ? R*(C+P)*S*H/(I*N2) : (dim==2) ? R*C*(S+P)*H/(I*N2) : R*C*S*(H+P)/(I*N2);
    const float z = 0.0f;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if (P==0) { cblas_scopy(R*C*S*H,X,1,Y,1); }
    else
    {
        for (int l=0, n=0; l<L; l++, n+=J)
        {
            cblas_scopy((N-P)*I,&X[n],1,&Y[n+P*I],1);
            cblas_scopy(P*I,&z,0,&Y[n],1);
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int prepad_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    if (R<0) { fprintf(stderr,"error in prepad_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in prepad_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in prepad_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in prepad_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * I;
    const int L = R*C*S*H/(I*N);
    const double z = 0.0;

    if (P==0) { cblas_dcopy(R*C*S*H,X,1,Y,1); }
    else if (P<=-N || P>=N) { cblas_dcopy(R*C*S*H,&z,0,Y,1); }
    else if (P<0)
    {
        for (int l=0, n=0; l<L; l++, n+=J)
        {
            cblas_dcopy((N+P)*I,&X[n-P*I],1,&Y[n],1);
            cblas_dcopy(-P*I,&z,0,&Y[n+(N+P)*I],1);
        }
    }
    else
    {
        for (int l=0, n=0; l<L; l++, n+=J)
        {
            cblas_dcopy((N-P)*I,&X[n],1,&Y[n+P*I],1);
            cblas_dcopy(P*I,&z,0,&Y[n],1);
        }
    }
    
    return 0;
}


int prepad_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    if (R<0) { fprintf(stderr,"error in prepad_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in prepad_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in prepad_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in prepad_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * I;
    const int L = R*C*S*H/(I*N);
    const float z[2] = {0.0f,0.0f};

    if (P==0) { cblas_ccopy(R*C*S*H,X,1,Y,1); }
    else if (P<=-N || P>=N) { cblas_ccopy(R*C*S*H,z,0,Y,1); }
    else if (P<0)
    {
        for (int l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((N+P)*I,&X[n-2*P*I],1,&Y[n],1);
            cblas_ccopy(-P*I,z,0,&Y[n+2*(N+P)*I],1);
        }
    }
    else
    {
        for (int l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((N-P)*I,&X[n],1,&Y[n+2*P*I],1);
            cblas_ccopy(P*I,z,0,&Y[n],1);
        }
    }
    
    return 0;
}


int prepad_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const int P)
{
    if (R<0) { fprintf(stderr,"error in prepad_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in prepad_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in prepad_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in prepad_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = N * I;
    const int L = R*C*S*H/(I*N);
    const double z[2] = {0.0,0.0};

    if (P==0) { cblas_zcopy(R*C*S*H,X,1,Y,1); }
    else if (P<=-N || P>=N) { cblas_zcopy(R*C*S*H,z,0,Y,1); }
    else if (P<0)
    {
        for (int l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((N+P)*I,&X[n-2*P*I],1,&Y[n],1);
            cblas_zcopy(-P*I,z,0,&Y[n+2*(N+P)*I],1);
        }
    }
    else
    {
        for (int l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((N-P)*I,&X[n],1,&Y[n+2*P*I],1);
            cblas_zcopy(P*I,z,0,&Y[n],1);
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
