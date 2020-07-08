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

int asum_s (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int asum_d (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int asum_c (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int asum_z (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);

int asum_inplace_s (float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int asum_inplace_d (double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int asum_inplace_c (float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int asum_inplace_z (double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);


int asum_s (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (size_t l=0; l<L; l++, X+=M*(N1-J))
    {
        for (size_t m=0; m<M; m++, X+=J, Y++)
        {
            *Y = cblas_sasum((int)N1,X,(int)K);
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int asum_d (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (size_t l=0; l<L; l++, X+=M*(N1-J))
    {
        for (size_t m=0; m<M; m++, X+=J, Y++)
        {
            *Y = cblas_dasum((int)N1,X,(int)K);
        }
    }

    return 0;
}


int asum_c (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (size_t l=0; l<L; l++, X+=2*M*(N1-J))
    {
        for (size_t m=0; m<M; m++, X+=2*J, Y++)
        {
            *Y = cblas_scasum((int)N1,X,(int)K);
        }
    }

    return 0;
}


int asum_z (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (size_t l=0; l<L; l++, X+=2*M*(N1-J))
    {
        for (size_t m=0; m<M; m++, X+=2*J, Y++)
        {
            *Y = cblas_dzasum((int)N1,X,(int)K);
        }
    }

    return 0;
}


int asum_inplace_s (float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (size_t l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
    {
        for (size_t m=0; m<M; m++, n+=J, n2++)
        {
            X[n2] = cblas_sasum((int)N1,&X[n],(int)K);
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int asum_inplace_d (double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (size_t l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
    {
        for (size_t m=0; m<M; m++, n+=J, n2++)
        {
            X[n2] = cblas_dasum((int)N1,&X[n],(int)K);
        }
    }

    return 0;
}


int asum_inplace_c (float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (size_t l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
    {
        for (size_t m=0; m<M; m++, n+=2*J, n2++)
        {
            X[n2] = cblas_scasum((int)N1,&X[n],(int)K);
        }
    }

    return 0;
}


int asum_inplace_z (double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    for (size_t l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
    {
        for (size_t m=0; m<M; m++, n+=2*J, n2++)
        {
            X[n2] = cblas_dzasum((int)N1,&X[n],(int)K);
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
