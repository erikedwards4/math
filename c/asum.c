//Gets sum of absolute-values (L1-norm) for each row or col of X according to dim.
//For complex case, output is real and is sum(|Xr|+|Xi|), not sum(|X|).
//The in-place version replaces the first N2 elements of X with the asum.

#include <stdio.h>
#include <cblas.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int asum_s (float *Y, const float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor);
int asum_d (double *Y, const double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor);
int asum_c (float *Y, const float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor);
int asum_z (double *Y, const double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor);

int asum_inplace_s (float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor);
int asum_inplace_d (double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor);
int asum_inplace_c (float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor);
int asum_inplace_z (double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor);


int asum_s (float *Y, const float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in asum_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in asum_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in asum_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in asum_s: H (num hyperslices X) must be positive\n"); return 1; }

    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = N/(M*N1);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (int l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
    {
        for (int m=0; m<M; m++, n+=J, n2++)
        {
            Y[n2] = cblas_sasum(N1,&X[n],K);
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int asum_d (double *Y, const double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in asum_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in asum_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in asum_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in asum_d: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = N/(M*N1);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (int l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
    {
        for (int m=0; m<M; m++, n+=J, n2++)
        {
            Y[n2] = cblas_dasum(N1,&X[n],K);
        }
    }

    return 0;
}


int asum_c (float *Y, const float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in asum_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in asum_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in asum_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in asum_c: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = N/(M*N1);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (int l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
    {
        for (int m=0; m<M; m++, n+=2*J, n2++)
        {
            Y[n2] = cblas_scasum(N1,&X[n],K);
        }
    }

    return 0;
}


int asum_z (double *Y, const double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in asum_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in asum_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in asum_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in asum_z: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = N/(M*N1);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (int l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
    {
        for (int m=0; m<M; m++, n+=2*J, n2++)
        {
            Y[n2] = cblas_dzasum(N1,&X[n],K);
        }
    }

    return 0;
}


int asum_inplace_s (float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in asum_inplace_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in asum_inplace_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in asum_inplace_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in asum_inplace_s: H (num hyperslices X) must be positive\n"); return 1; }

    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = N/(M*N1);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (int l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
    {
        for (int m=0; m<M; m++, n+=J, n2++)
        {
            X[n2] = cblas_sasum(N1,&X[n],K);
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int asum_inplace_d (double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in asum_inplace_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in asum_inplace_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in asum_inplace_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in asum_inplace_d: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = N/(M*N1);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (int l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
    {
        for (int m=0; m<M; m++, n+=J, n2++)
        {
            X[n2] = cblas_dasum(N1,&X[n],K);
        }
    }

    return 0;
}


int asum_inplace_c (float *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in asum_inplace_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in asum_inplace_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in asum_inplace_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in asum_inplace_c: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = N/(M*N1);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (int l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
    {
        for (int m=0; m<M; m++, n+=2*J, n2++)
        {
            X[n2] = cblas_scasum(N1,&X[n],K);
        }
    }

    return 0;
}


int asum_inplace_z (double *X, const int R, const int C,const int S, const int H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in asum_inplace_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in asum_inplace_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in asum_inplace_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in asum_inplace_z: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = N/(M*N1);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (int l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
    {
        for (int m=0; m<M; m++, n+=2*J, n2++)
        {
            X[n2] = cblas_dzasum(N1,&X[n],K);
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
