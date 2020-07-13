//Gets count of nonzero values for each row or col of X according to dim.
//This is also the Hamming norm (or L0 norm) of each vector in X.

//For complex case, real and imag parts calculated separately.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cnt_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor);
int cnt_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor);
int cnt_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor);
int cnt_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor);

int cnt_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor);
int cnt_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor);
int cnt_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor);
int cnt_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor);


int cnt_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
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
            *Y = 0.0f;
            for (size_t n1=0; n1<N1; n1++) { *Y += (float)(*(X+n1*K)!=0.0f); }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int cnt_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
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
            *Y = 0.0;
            for (size_t n1=0; n1<N1; n1++) { *Y += (double)(*(X+n1*K)!=0.0); }
        }
    }

    return 0;
}


int cnt_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
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
            *Y = 0.0f;
            for (size_t n1=0; n1<2*N1; n1+=2) { *Y += (float)(*(X+n1*K)!=0.0f || *(X+n1*K+1)!=0.0f); }
        }
    }

    return 0;
}


int cnt_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
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
            *Y = 0.0;
            for (size_t n1=0; n1<2*N1; n1+=2) { *Y += (double)(*(X+n1*K)!=0.0 || *(X+n1*K+1)!=0.0); }
        }
    }

    return 0;
}


int cnt_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
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
            float sm = 0.0f;
            for (size_t n1=0; n1<N1; n1++) { sm += (float)(X[n+n1*K]!=0.0f); }
            X[n2] = sm;
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int cnt_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
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
            double sm = 0.0;
            for (size_t n1=0; n1<N1; n1++) { sm += (double)(X[n+n1*K]!=0.0); }
            X[n2] = sm;
        }
    }

    return 0;
}


int cnt_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
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
            float sm = 0.0f;
            for (size_t n1=0; n1<2*N1; n1+=2) { sm += (float)(X[n+n1*K]!=0.0f || X[n+n1*K+1]!=0.0f); }
            X[n2] = sm;
        }
    }

    return 0;
}


int cnt_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
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
            double sm = 0.0;
            for (size_t n1=0; n1<2*N1; n1+=2) { sm += (double)(X[n+n1*K]!=0.0 || X[n+n1*K+1]!=0.0); }
            X[n2] = sm;
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
