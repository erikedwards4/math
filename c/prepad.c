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

int prepad_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P);
int prepad_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P);
int prepad_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P);
int prepad_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P);


int prepad_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P)
{
    if (dim>3) { fprintf(stderr,"error in prepad_s: dim must be in [0 3]\n"); return 1; }

    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t N2 = N1 + P;
    const size_t I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J1 = N1*I, J2 = N2*I;
    const size_t L1 = R*C*S*H/(I*N1);
    const size_t L2 = (dim==0) ? (R+P)*C*S*H/(I*N2) : (dim==1) ? R*(C+P)*S*H/(I*N2) : (dim==2) ? R*C*(S+P)*H/(I*N2) : R*C*S*(H+P)/(I*N2);
    const float z = 0.0f;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if (P==0) { cblas_scopy((int)(R*C*S*H),X,1,Y,1); }
    else
    {
        for (size_t l=0, n=0; l<L; ++l, n+=J)
        {
            cblas_scopy((int)((N-P)*I),&X[n],1,&Y[n+P*I],1);
            cblas_scopy((int)(P*I)&z,0,&Y[n],1);
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int prepad_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P)
{
    if (dim>3) { fprintf(stderr,"error in prepad_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N * I;
    const size_t L = R*C*S*H/(I*N);
    const double z = 0.0;

    if (P==0) { cblas_dcopy((int)(R*C*S*H),X,1,Y,1); }
    else if (P<=-N || P>=N) { cblas_dcopy((int)(R*C*S*H),&z,0,Y,1); }
    else if (P<0)
    {
        for (size_t l=0, n=0; l<L; ++l, n+=J)
        {
            cblas_dcopy((N+P)*I,&X[n-P*I],1,&Y[n],1);
            cblas_dcopy((int)(-P*I)&z,0,&Y[n+(N+P)*I],1);
        }
    }
    else
    {
        for (size_t l=0, n=0; l<L; ++l, n+=J)
        {
            cblas_dcopy((int)((N-P)*I),&X[n],1,&Y[n+P*I],1);
            cblas_dcopy((int)(P*I)&z,0,&Y[n],1);
        }
    }
    
    return 0;
}


int prepad_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P)
{
    if (dim>3) { fprintf(stderr,"error in prepad_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N * I;
    const size_t L = R*C*S*H/(I*N);
    const float z[2] = {0.0f,0.0f};

    if (P==0) { cblas_ccopy((int)(R*C*S*H),X,1,Y,1); }
    else if (P<=-N || P>=N) { cblas_ccopy((int)(R*C*S*H),z,0,Y,1); }
    else if (P<0)
    {
        for (size_t l=0, n=0; l<L; ++l, n+=2*J)
        {
            cblas_ccopy((N+P)*I,&X[n-2*P*I],1,&Y[n],1);
            cblas_ccopy((int)(-P*I)z,0,&Y[n+2*(N+P)*I],1);
        }
    }
    else
    {
        for (size_t l=0, n=0; l<L; ++l, n+=2*J)
        {
            cblas_ccopy((int)((N-P)*I),&X[n],1,&Y[n+2*P*I],1);
            cblas_ccopy((int)(P*I)z,0,&Y[n],1);
        }
    }
    
    return 0;
}


int prepad_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P)
{
    if (dim>3) { fprintf(stderr,"error in prepad_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N * I;
    const size_t L = R*C*S*H/(I*N);
    const double z[2] = {0.0,0.0};

    if (P==0) { cblas_zcopy((int)(R*C*S*H),X,1,Y,1); }
    else if (P<=-N || P>=N) { cblas_zcopy((int)(R*C*S*H),z,0,Y,1); }
    else if (P<0)
    {
        for (size_t l=0, n=0; l<L; ++l, n+=2*J)
        {
            cblas_zcopy((N+P)*I,&X[n-2*P*I],1,&Y[n],1);
            cblas_zcopy((int)(-P*I)z,0,&Y[n+2*(N+P)*I],1);
        }
    }
    else
    {
        for (size_t l=0, n=0; l<L; ++l, n+=2*J)
        {
            cblas_zcopy((int)((N-P)*I),&X[n],1,&Y[n+2*P*I],1);
            cblas_zcopy((int)(P*I)z,0,&Y[n],1);
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
