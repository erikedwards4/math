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

int cshift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P);
int cshift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P);
int cshift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P);
int cshift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P);

int cshift_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P);
int cshift_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P);
int cshift_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P);
int cshift_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P);


int cshift_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P)
{
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);
    const int D = P % (int)N1;

    if (D==0) { cblas_scopy((int)(R*C*S*H),X,1,Y,1); }
    else if (D<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_scopy((int)K*((int)N1+D),&X[(size_t)((int)n-D*(int)K)],1,&Y[n],1);
            cblas_scopy(-D*(int)K,&X[n],1,&Y[n+(size_t)((int)N1+D)*K],1);
        }
    }
    else
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_scopy((int)K*((int)N1-D),&X[n],1,&Y[n+(size_t)D*K],1);
            cblas_scopy(D*(int)K,&X[n+(N1-(size_t)D)*K],1,&Y[n],1);
        }
    }

    return 0;
}


int cshift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P)
{
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);
    const int D = P % (int)N1;

    if (D==0) { cblas_dcopy((int)(R*C*S*H),X,1,Y,1); }
    else if (D<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_dcopy((int)K*((int)N1+D),&X[(size_t)((int)n-D*(int)K)],1,&Y[n],1);
            cblas_dcopy(-D*(int)K,&X[n],1,&Y[n+(size_t)((int)N1+D)*K],1);
        }
    }
    else
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_dcopy((int)K*((int)N1-D),&X[n],1,&Y[n+(size_t)D*K],1);
            cblas_dcopy(D*(int)K,&X[n+(N1-(size_t)D)*K],1,&Y[n],1);
        }
    }

    return 0;
}


int cshift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P)
{
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);
    const int D = P % (int)N1;

    if (D==0) { cblas_ccopy((int)(R*C*S*H),X,1,Y,1); }
    else if (D<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((int)K*((int)N1+D),&X[(size_t)((int)n-2*D*(int)K)],1,&Y[n],1);
            cblas_ccopy(-D*(int)K,&X[n],1,&Y[n+2*(size_t)((int)N1+D)*K],1);
        }
    }
    else
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((int)K*((int)N1-D),&X[n],1,&Y[n+2*(size_t)D*K],1);
            cblas_ccopy(D*(int)K,&X[n+2*(N1-(size_t)D)*K],1,&Y[n],1);
        }
    }

    return 0;
}


int cshift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P)
{
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);
    const int D = P % (int)N1;

    if (D==0) { cblas_zcopy((int)(R*C*S*H),X,1,Y,1); }
    else if (D<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((int)K*((int)N1+D),&X[(size_t)((int)n-2*D*(int)K)],1,&Y[n],1);
            cblas_zcopy(-D*(int)K,&X[n],1,&Y[n+2*(size_t)((int)N1+D)*K],1);
        }
    }
    else
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((int)K*((int)N1-D),&X[n],1,&Y[n+2*(size_t)D*K],1);
            cblas_zcopy(D*(int)K,&X[n+2*(N1-(size_t)D)*K],1,&Y[n],1);
        }
    }

    return 0;
}


int cshift_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P)
{
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);
    const int D = P % (int)N1;
    const size_t MS = (D>0) ? K*(N1-(size_t)D) : (size_t)(-D*(int)K);

    float *X1;
    if (!(X1=(float *)malloc(MS*sizeof(float)))) { fprintf(stderr,"error in cshift_inplace_s: problem with malloc. "); perror("malloc"); return 1; }

    if (D<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_scopy(-D*(int)K,&X[n],1,X1,1);
            cblas_scopy((int)K*((int)N1+D),&X[(size_t)((int)n-D*(int)K)],1,&X[n],1);
            cblas_scopy(-D*(int)K,X1,1,&X[n+(size_t)((int)N1+D)*K],1);
        }
    }
    else if (D>0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_scopy((int)K*((int)N1-D),&X[n],1,X1,1);
            cblas_scopy(D*(int)K,&X[n+(N1-(size_t)D)*K],1,&X[n],1);
            cblas_scopy((int)K*((int)N1-D),X1,1,&X[n+(size_t)D*K],1);
        }
    }

    free(X1);
    return 0;
}


int cshift_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P)
{
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);
    const int D = P % (int)N1;
    const size_t MS = (D>0) ? K*(N1-(size_t)D) : (size_t)(-D*(int)K);

    double *X1;
    if (!(X1=(double *)malloc(MS*sizeof(double)))) { fprintf(stderr,"error in cshift_inplace_d: problem with malloc. "); perror("malloc"); return 1; }

    if (D<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_dcopy(-D*(int)K,&X[n],1,X1,1);
            cblas_dcopy((int)K*((int)N1+D),&X[(size_t)((int)n-D*(int)K)],1,&X[n],1);
            cblas_dcopy(-D*(int)K,X1,1,&X[n+(size_t)((int)N1+D)*K],1);
        }
    }
    else if (D>0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=J)
        {
            cblas_dcopy((int)K*((int)N1-D),&X[n],1,X1,1);
            cblas_dcopy(D*(int)K,&X[n+(N1-(size_t)D)*K],1,&X[n],1);
            cblas_dcopy((int)K*((int)N1-D),X1,1,&X[n+(size_t)D*K],1);
        }
    }

    free(X1);
    return 0;
}


int cshift_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P)
{
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);
    const int D = P % (int)N1;
    const size_t MS = (D>0) ? K*(N1-(size_t)D) : (size_t)(-D*(int)K);

    float *X1;
    if (!(X1=(float *)malloc(MS*2u*sizeof(float)))) { fprintf(stderr,"error in cshift_inplace_c: problem with malloc. "); perror("malloc"); return 1; }

    if (D<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy(-D*(int)K,&X[n],1,X1,1);
            cblas_ccopy((int)K*((int)N1+D),&X[(size_t)((int)n-2*D*(int)K)],1,&X[n],1);
            cblas_ccopy(-D*(int)K,X1,1,&X[n+2*(size_t)((int)N1+D)*K],1);
        }
    }
    else if (D>0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_ccopy((int)K*((int)N1-D),&X[n],1,X1,1);
            cblas_ccopy(D*(int)K,&X[n+2*(N1-(size_t)D)*K],1,&X[n],1);
            cblas_ccopy((int)K*((int)N1-D),X1,1,&X[n+2*(size_t)D*K],1);
        }
    }

    free(X1);
    return 0;
}


int cshift_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int dim, const int P)
{
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = N1 * K;
    const size_t L = R*C*S*H/(K*N1);
    const int D = P % (int)N1;
    const size_t MS = (D>0) ? K*(N1-(size_t)D) : (size_t)(-D*(int)K);

    double *X1;
    if (!(X1=(double *)malloc(MS*2u*sizeof(double)))) { fprintf(stderr,"error in cshift_inplace_z: problem with malloc. "); perror("malloc"); return 1; }

    if (D<0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy(-D*(int)K,&X[n],1,X1,1);
            cblas_zcopy((int)K*((int)N1+D),&X[(size_t)((int)n-2*D*(int)K)],1,&X[n],1);
            cblas_zcopy(-D*(int)K,X1,1,&X[n+2*(size_t)((int)N1+D)*K],1);
        }
    }
    else if (D>0)
    {
        for (size_t l=0, n=0; l<L; l++, n+=2*J)
        {
            cblas_zcopy((int)K*((int)N1-D),&X[n],1,X1,1);
            cblas_zcopy(D*(int)K,&X[n+2*(N1-(size_t)D)*K],1,&X[n],1);
            cblas_zcopy((int)K*((int)N1-D),X1,1,&X[n+2*(size_t)D*K],1);
        }
    }

    free(X1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
