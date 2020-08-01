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
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t L1 = (P<0) ? K*(size_t)((int)L+P) : K*(L-(size_t)P);
    const size_t L2 = (P<0) ? K*(size_t)(-P) : K*(size_t)P;

    if (N==0) {}
    else if (P==0) { cblas_scopy((int)N,X,1,Y,1); }
    else if (P<=-(int)L || P>=(int)L) { cblas_scopy((int)N,&z,0,Y,1); }
    else if (L==N)
    {
        if (P<0)
        {
            cblas_scopy((int)L1,&X[L2],1,Y,1);
            cblas_scopy((int)L2,&z,0,&Y[L1],1);
        }
        else
        {
            cblas_scopy((int)L2,&z,0,Y,1);
            cblas_scopy((int)L1,X,1,&Y[L2],1);
        }
    }
    else
    {
        const size_t J = L*K, G = N/J;
        if (P<0)
        {
            for (size_t g=0, n=0; g<G; ++g, n+=J)
            {
                cblas_scopy((int)L1,&X[n+L2],1,&Y[n],1);
                cblas_scopy((int)L2,&z,0,&Y[n+L1],1);
            }
        }
        else
        {
            for (size_t g=0, n=0; g<G; ++g, n+=J)
            {
                cblas_scopy((int)L2,&z,0,&Y[n],1);
                cblas_scopy((int)L1,&X[n],1,&Y[n+L2],1);
            }
        }
    }
    
    return 0;
}


int shift_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_d: dim must be in [0 3]\n"); return 1; }

    const double z = 0.0;
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t L1 = (P<0) ? K*(size_t)((int)L+P) : K*(L-(size_t)P);
    const size_t L2 = (P<0) ? K*(size_t)(-P) : K*(size_t)P;

    if (N==0) {}
    else if (P==0) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (P<=-(int)L || P>=(int)L) { cblas_dcopy((int)N,&z,0,Y,1); }
    else if (L==N)
    {
        if (P<0)
        {
            cblas_dcopy((int)L1,&X[L2],1,Y,1);
            cblas_dcopy((int)L2,&z,0,&Y[L1],1);
        }
        else
        {
            cblas_dcopy((int)L2,&z,0,Y,1);
            cblas_dcopy((int)L1,X,1,&Y[L2],1);
        }
    }
    else
    {
        const size_t J = L*K, G = N/J;
        if (P<0)
        {
            for (size_t g=0, n=0; g<G; ++g, n+=J)
            {
                cblas_dcopy((int)L1,&X[n+L2],1,&Y[n],1);
                cblas_dcopy((int)L2,&z,0,&Y[n+L1],1);
            }
        }
        else
        {
            for (size_t g=0, n=0; g<G; ++g, n+=J)
            {
                cblas_dcopy((int)L2,&z,0,&Y[n],1);
                cblas_dcopy((int)L1,&X[n],1,&Y[n+L2],1);
            }
        }
    }
    
    return 0;
}


int shift_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_c: dim must be in [0 3]\n"); return 1; }

    const float z[2] = {0.0f,0.0f};
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t L1 = (P<0) ? K*(size_t)((int)L+P) : K*(L-(size_t)P);
    const size_t L2 = (P<0) ? K*(size_t)(-P) : K*(size_t)P;

    if (N==0) {}
    else if (P==0) { cblas_ccopy((int)N,X,1,Y,1); }
    else if (P<=-(int)L || P>=(int)L) { cblas_ccopy((int)N,z,0,Y,1); }
    else if (L==N)
    {
        if (P<0)
        {
            cblas_ccopy((int)L1,&X[L2],1,Y,1);
            cblas_ccopy((int)L2,z,0,&Y[L1],1);
        }
        else
        {
            cblas_ccopy((int)L2,z,0,Y,1);
            cblas_ccopy((int)L1,X,1,&Y[L2],1);
        }
    }
    else
    {
        const size_t J = L*K, G = N/J;
        if (P<0)
        {
            for (size_t g=0, n=0; g<G; ++g, n+=2*J)
            {
                cblas_ccopy((int)L1,&X[n+2*L2],1,&Y[n],1);
                cblas_ccopy((int)L2,z,0,&Y[n+2*L1],1);
            }
        }
        else
        {
            for (size_t g=0, n=0; g<G; ++g, n+=2*J)
            {
                cblas_ccopy((int)L2,z,0,&Y[n],1);
                cblas_ccopy((int)L1,&X[n],1,&Y[n+2*L2],1);
            }
        }
    }
    
    return 0;
}


int shift_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_z: dim must be in [0 3]\n"); return 1; }

    const double z[2] = {0.0,0.0};
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t L1 = (P<0) ? K*(size_t)((int)L+P) : K*(L-(size_t)P);
    const size_t L2 = (P<0) ? K*(size_t)(-P) : K*(size_t)P;

    if (N==0) {}
    else if (P==0) { cblas_zcopy((int)N,X,1,Y,1); }
    else if (P<=-(int)L || P>=(int)L) { cblas_zcopy((int)N,z,0,Y,1); }
    else if (L==N)
    {
        if (P<0)
        {
            cblas_zcopy((int)L1,&X[L2],1,Y,1);
            cblas_zcopy((int)L2,z,0,&Y[L1],1);
        }
        else
        {
            cblas_zcopy((int)L2,z,0,Y,1);
            cblas_zcopy((int)L1,X,1,&Y[L2],1);
        }
    }
    else
    {
        const size_t J = L*K, G = N/J;
        if (P<0)
        {
            for (size_t g=0, n=0; g<G; ++g, n+=2*J)
            {
                cblas_zcopy((int)L1,&X[n+2*L2],1,&Y[n],1);
                cblas_zcopy((int)L2,z,0,&Y[n+2*L1],1);
            }
        }
        else
        {
            for (size_t g=0, n=0; g<G; ++g, n+=2*J)
            {
                cblas_zcopy((int)L2,z,0,&Y[n],1);
                cblas_zcopy((int)L1,&X[n],1,&Y[n+2*L2],1);
            }
        }
    }
    
    return 0;
}


int shift_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_inplace_s: dim must be in [0 3]\n"); return 1; }

    const float z = 0.0f;
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    if (N==0 || P==0) {}
    else if (P<=-(int)L || P>=(int)L) { cblas_scopy((int)N,&z,0,X,1); }
    else if (P<0)
    {
        const size_t L1 = (P<0) ? K*(size_t)((int)L+P) : K*(L-(size_t)P);
        const size_t L2 = (P<0) ? K*(size_t)(-P) : K*(size_t)P;
        if (L==N)
        {
            cblas_scopy((int)L1,&X[L2],1,X,1);
            cblas_scopy((int)L2,&z,0,&X[L1],1);
        }
        else
        {
            const size_t J = L*K, G = N/J;
            for (size_t g=0, n=0; g<G; ++g, n+=J)
            {
                cblas_scopy((int)L1,&X[n+L2],1,&X[n],1);
                cblas_scopy((int)L2,&z,0,&X[n+L1],1);
            }
        }
    }
    else
    {
        const size_t L1 = (P<0) ? K*(size_t)((int)L+P) : K*(L-(size_t)P);
        const size_t L2 = (P<0) ? K*(size_t)(-P) : K*(size_t)P;
        const size_t MS = (P<0) ? L2 : L1;
        float *X1;
        if (!(X1=(float *)malloc(MS*sizeof(float)))) { fprintf(stderr,"error in shift_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        if (L==N)
        {
            cblas_scopy((int)L1,X,1,X1,1);
            cblas_scopy((int)L2,&z,0,X,1);
            cblas_scopy((int)L1,X1,1,&X[L2],1);
        }
        else
        {
            const size_t J = L*K, G = N/J;
            for (size_t g=0, n=0; g<G; ++g, n+=J)
            {
                cblas_scopy((int)L1,&X[n],1,X1,1);
                cblas_scopy((int)L2,&z,0,&X[n],1);
                cblas_scopy((int)L1,X1,1,&X[n+L2],1);
            }
        }
        free(X1);
    }
    
    return 0;
}


int shift_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_inplace_d: dim must be in [0 3]\n"); return 1; }

    const double z = 0.0;
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    if (N==0 || P==0) {}
    else if (P<=-(int)L || P>=(int)L) { cblas_dcopy((int)N,&z,0,X,1); }
    else if (P<0)
    {
        const size_t L1 = (P<0) ? K*(size_t)((int)L+P) : K*(L-(size_t)P);
        const size_t L2 = (P<0) ? K*(size_t)(-P) : K*(size_t)P;
        if (L==N)
        {
            cblas_dcopy((int)L1,&X[L2],1,X,1);
            cblas_dcopy((int)L2,&z,0,&X[L1],1);
        }
        else
        {
            const size_t J = L*K, G = N/J;
            for (size_t g=0, n=0; g<G; ++g, n+=J)
            {
                cblas_dcopy((int)L1,&X[n+L2],1,&X[n],1);
                cblas_dcopy((int)L2,&z,0,&X[n+L1],1);
            }
        }
    }
    else
    {
        const size_t L1 = (P<0) ? K*(size_t)((int)L+P) : K*(L-(size_t)P);
        const size_t L2 = (P<0) ? K*(size_t)(-P) : K*(size_t)P;
        const size_t MS = (P<0) ? L2 : L1;
        double *X1;
        if (!(X1=(double *)malloc(MS*sizeof(double)))) { fprintf(stderr,"error in shift_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        if (L==N)
        {
            cblas_dcopy((int)L1,X,1,X1,1);
            cblas_dcopy((int)L2,&z,0,X,1);
            cblas_dcopy((int)L1,X1,1,&X[L2],1);
        }
        else
        {
            const size_t J = L*K, G = N/J;
            for (size_t g=0, n=0; g<G; ++g, n+=J)
            {
                cblas_dcopy((int)L1,&X[n],1,X1,1);
                cblas_dcopy((int)L2,&z,0,&X[n],1);
                cblas_dcopy((int)L1,X1,1,&X[n+L2],1);
            }
        }
        free(X1);
    }
    
    return 0;
}


int shift_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_inplace_c: dim must be in [0 3]\n"); return 1; }

    const float z[2] = {0.0f,0.0f};
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    if (N==0 || P==0) {}
    else if (P<=-(int)L || P>=(int)L) { cblas_ccopy((int)N,z,0,X,1); }
    else if (P<0)
    {
        const size_t L1 = (P<0) ? K*(size_t)((int)L+P) : K*(L-(size_t)P);
        const size_t L2 = (P<0) ? K*(size_t)(-P) : K*(size_t)P;
        if (L==N)
        {
            cblas_ccopy((int)L1,&X[2*L2],1,X,1);
            cblas_ccopy((int)L2,z,0,&X[2*L1],1);
        }
        else
        {
            const size_t J = L*K, G = N/J;
            for (size_t g=0, n=0; g<G; ++g, n+=2*J)
            {
                cblas_ccopy((int)L1,&X[n+2*L2],1,&X[n],1);
                cblas_ccopy((int)L2,z,0,&X[n+2*L1],1);
            }
        }
    }
    else
    {
        const size_t L1 = (P<0) ? K*(size_t)((int)L+P) : K*(L-(size_t)P);
        const size_t L2 = (P<0) ? K*(size_t)(-P) : K*(size_t)P;
        const size_t MS = (P<0) ? L2 : L1;
        float *X1;
        if (!(X1=(float *)malloc(2*MS*sizeof(float)))) { fprintf(stderr,"error in shift_inplace_c: problem with malloc. "); perror("malloc"); return 1; }
        if (L==N)
        {
            cblas_ccopy((int)L1,X,1,X1,1);
            cblas_ccopy((int)L2,z,0,X,1);
            cblas_ccopy((int)L1,X1,1,&X[2*L2],1);
        }
        else
        {
            const size_t J = L*K, G = N/J;
            for (size_t g=0, n=0; g<G; ++g, n+=2*J)
            {
                cblas_ccopy((int)L1,&X[n],1,X1,1);
                cblas_ccopy((int)L2,z,0,&X[n],1);
                cblas_ccopy((int)L1,X1,1,&X[n+2*L2],1);
            }
        }
        free(X1);
    }
    
    return 0;
}


int shift_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int P)
{
    if (dim>3) { fprintf(stderr,"error in shift_inplace_z: dim must be in [0 3]\n"); return 1; }

    const double z[2] = {0.0,0.0};
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    if (N==0 || P==0) {}
    else if (P<=-(int)L || P>=(int)L) { cblas_zcopy((int)N,z,0,X,1); }
    else if (P<0)
    {
        const size_t L1 = (P<0) ? K*(size_t)((int)L+P) : K*(L-(size_t)P);
        const size_t L2 = (P<0) ? K*(size_t)(-P) : K*(size_t)P;
        if (L==N)
        {
            cblas_zcopy((int)L1,&X[2*L2],1,X,1);
            cblas_zcopy((int)L2,z,0,&X[2*L1],1);
        }
        else
        {
            const size_t J = L*K, G = N/J;
            for (size_t g=0, n=0; g<G; ++g, n+=2*J)
            {
                cblas_zcopy((int)L1,&X[n+2*L2],1,&X[n],1);
                cblas_zcopy((int)L2,z,0,&X[n+2*L1],1);
            }
        }
    }
    else
    {
        const size_t L1 = (P<0) ? K*(size_t)((int)L+P) : K*(L-(size_t)P);
        const size_t L2 = (P<0) ? K*(size_t)(-P) : K*(size_t)P;
        const size_t MS = (P<0) ? L2 : L1;
        double *X1;
        if (!(X1=(double *)malloc(2*MS*sizeof(double)))) { fprintf(stderr,"error in shift_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
        if (L==N)
        {
            cblas_zcopy((int)L1,X,1,X1,1);
            cblas_zcopy((int)L2,z,0,X,1);
            cblas_zcopy((int)L1,X1,1,&X[2*L2],1);
        }
        else
        {
            const size_t J = L*K, G = N/J;
            for (size_t g=0, n=0; g<G; ++g, n+=2*J)
            {
                cblas_zcopy((int)L1,&X[n],1,X1,1);
                cblas_zcopy((int)L2,z,0,&X[n],1);
                cblas_zcopy((int)L1,X1,1,&X[n+2*L2],1);
            }
        }
        free(X1);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
