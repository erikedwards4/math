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

int shift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int shift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int shift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int shift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);

int shift_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int shift_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int shift_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);
int shift_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P);


int shift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_s: dim must be in [0 3]\n"); return 1; }

    const float z = 0.0f;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);

    if (P==0) { cblas_scopy((int)(R*C*S*H),X,1,Y,1); }
    else if (P<=-(int)N1 || P>=(int)N1) { cblas_scopy((int)(R*C*S*H),&z,0,Y,1); }
    else if (P<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_scopy((int)K*((int)N1+P),&X[(size_t)((int)n-P*(int)K)],1,&Y[n],1);
            cblas_scopy(-P*(int)K,&z,0,&Y[n+(size_t)((int)N1+P)*K],1);
        }
    }
    else
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_scopy((int)K*((int)N1-P),&X[n],1,&Y[n+(size_t)P*K],1);
            cblas_scopy(P*(int)K,&z,0,&Y[n],1);
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int shift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_d: dim must be in [0 3]\n"); return 1; }

    const double z = 0.0;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);

    if (P==0) { cblas_dcopy((int)(R*C*S*H),X,1,Y,1); }
    else if (P<=-(int)N1 || P>=(int)N1) { cblas_dcopy((int)(R*C*S*H),&z,0,Y,1); }
    else if (P<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_dcopy((int)K*((int)N1+P),&X[(size_t)((int)n-P*(int)K)],1,&Y[n],1);
            cblas_dcopy(-P*(int)K,&z,0,&Y[n+(size_t)((int)N1+P)*K],1);
        }
    }
    else
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_dcopy((int)K*((int)N1-P),&X[n],1,&Y[n+(size_t)P*K],1);
            cblas_dcopy(P*(int)K,&z,0,&Y[n],1);
        }
    }
    
    return 0;
}


int shift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_c: dim must be in [0 3]\n"); return 1; }

    const float z[2] = {0.0f,0.0f};
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);

    if (P==0) { cblas_ccopy((int)(R*C*S*H),X,1,Y,1); }
    else if (P<=-(int)N1 || P>=(int)N1) { cblas_ccopy((int)(R*C*S*H),z,0,Y,1); }
    else if (P<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((int)K*((int)N1+P),&X[(size_t)((int)n-2*P*(int)K)],1,&Y[n],1);
            cblas_ccopy(-P*(int)K,z,0,&Y[n+2*(size_t)((int)N1+P)*K],1);
        }
    }
    else
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((int)K*((int)N1-P),&X[n],1,&Y[n+2*(size_t)P*K],1);
            cblas_ccopy(P*(int)K,z,0,&Y[n],1);
        }
    }
    
    return 0;
}


int shift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_z: dim must be in [0 3]\n"); return 1; }

    const double z[2] = {0.0,0.0};
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);

    if (P==0) { cblas_zcopy((int)(R*C*S*H),X,1,Y,1); }
    else if (P<=-(int)N1 || P>=(int)N1) { cblas_zcopy((int)(R*C*S*H),z,0,Y,1); }
    else if (P<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((int)K*((int)N1+P),&X[(size_t)((int)n-2*P*(int)K)],1,&Y[n],1);
            cblas_zcopy(-P*(int)K,z,0,&Y[n+2*(size_t)((int)N1+P)*K],1);
        }
    }
    else
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((int)K*((int)N1-P),&X[n],1,&Y[n+2*(size_t)P*K],1);
            cblas_zcopy(P*(int)K,z,0,&Y[n],1);
        }
    }
    
    return 0;
}


int shift_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_inplace_s: dim must be in [0 3]\n"); return 1; }

    const float z = 0.0f;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);

    if (P<=-(int)N1 || P>=(int)N1) { cblas_scopy((int)(R*C*S*H),&z,0,X,1); }
    else if (P<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_scopy((int)K*((int)N1+P),&X[(size_t)((int)n-P*(int)K)],1,&X[n],1);
            cblas_scopy(-P*(int)K,&z,0,&X[n+(size_t)((int)N1+P)*K],1);
        }
    }
    else if (P>0)
    {
        float *X1;
        if (!(X1=(float *)malloc(K*(N1-(size_t)P)*sizeof(float)))) { fprintf(stderr,"error in shift_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_scopy((int)K*((int)N1-P),&X[n],1,X1,1);
            cblas_scopy((int)K*((int)N1-P),X1,1,&X[n+(size_t)P*K],1);
            cblas_scopy(P*(int)K,&z,0,&X[n],1);
        }
        free(X1);
    }

    return 0;
}


int shift_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_inplace_d: dim must be in [0 3]\n"); return 1; }

    const double z = 0.0;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);

    if (P<=-(int)N1 || P>=(int)N1) { cblas_dcopy((int)(R*C*S*H),&z,0,X,1); }
    else if (P<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_dcopy((int)K*((int)N1+P),&X[(size_t)((int)n-P*(int)K)],1,&X[n],1);
            cblas_dcopy(-P*(int)K,&z,0,&X[n+(size_t)((int)N1+P)*K],1);
        }
    }
    else if (P>0)
    {
        double *X1;
        if (!(X1=(double *)malloc(K*(N1-(size_t)P)*sizeof(double)))) { fprintf(stderr,"error in shift_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_dcopy((int)K*((int)N1-P),&X[n],1,X1,1);
            cblas_dcopy((int)K*((int)N1-P),X1,1,&X[n+(size_t)P*K],1);
            cblas_dcopy(P*(int)K,&z,0,&X[n],1);
        }
        free(X1);
    }

    return 0;
}


int shift_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_inplace_c: dim must be in [0 3]\n"); return 1; }

    const float z[2] = {0.0f,0.0f};
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);

    if (P<=-(int)N1 || P>=(int)N1) { cblas_ccopy((int)(R*C*S*H),z,0,X,1); }
    else if (P<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((int)K*((int)N1+P),&X[(size_t)((int)n-2*P*(int)K)],1,&X[n],1);
            cblas_ccopy(-P*(int)K,z,0,&X[n+2*(size_t)((int)N1+P)*K],1);
        }
    }
    else if (P>0)
    {
        float *X1;
        if (!(X1=(float *)malloc(2*K*(N1-(size_t)P)*sizeof(float)))) { fprintf(stderr,"error in shift_inplace_c: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((int)K*((int)N1-P),&X[n],1,X1,1);
            cblas_ccopy((int)K*((int)N1-P),X1,1,&X[n+2*(size_t)P*K],1);
            cblas_ccopy(P*(int)K,z,0,&X[n],1);
        }
        free(X1);
    }
    
    return 0;
}


int shift_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_inplace_z: dim must be in [0 3]\n"); return 1; }

    const double z[2] = {0.0,0.0};
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);

    if (P<=-(int)N1 || P>=(int)N1) { cblas_zcopy((int)(R*C*S*H),z,0,X,1); }
    else if (P<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((int)K*((int)N1+P),&X[(size_t)((int)n-2*P*(int)K)],1,&X[n],1);
            cblas_zcopy(-P*(int)K,z,0,&X[n+2*(size_t)((int)N1+P)*K],1);
        }
    }
    else if (P>0)
    {
        double *X1;
        if (!(X1=(double *)malloc(2*K*(N1-(size_t)P)*sizeof(double)))) { fprintf(stderr,"error in shift_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((int)K*((int)N1-P),&X[n],1,X1,1);
            cblas_zcopy((int)K*((int)N1-P),X1,1,&X[n+2*(size_t)P*K],1);
            cblas_zcopy(P*(int)K,z,0,&X[n],1);
        }
        free(X1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
