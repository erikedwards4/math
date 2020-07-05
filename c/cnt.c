//Gets count of nonzero values for each row or col of X according to dim.
//For complex case, real and imag parts calculated separately.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cnt_s (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int cnt_d (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int cnt_c (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int cnt_z (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);

int cnt_inplace_s (float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int cnt_inplace_d (double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int cnt_inplace_c (float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int cnt_inplace_z (double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);


int cnt_s (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in cnt_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in cnt_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in cnt_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in cnt_s: H (num hyperslices X) must be positive\n"); return 1; }

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
            Y[n2] = 0.0f;
            for (int n1=0; n1<N1; n1++) { Y[n2] += (float)(X[n+n1*K]!=0.0f); }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int cnt_d (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in cnt_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in cnt_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in cnt_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in cnt_d: H (num hyperslices X) must be positive\n"); return 1; }

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
            Y[n2] = 0.0;
            for (int n1=0; n1<N1; n1++) { Y[n2] += (double)(X[n+n1*K]!=0.0); }
        }
    }

    return 0;
}


int cnt_c (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in cnt_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in cnt_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in cnt_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in cnt_c: H (num hyperslices X) must be positive\n"); return 1; }

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
            Y[n2] = 0.0f;
            for (int n1=0; n1<2*N1; n1+=2) { Y[n2] += (float)(X[n+n1*K]!=0.0f || X[n+n1*K+1]!=0.0f); }
        }
    }

    return 0;
}


int cnt_z (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in cnt_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in cnt_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in cnt_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in cnt_z: H (num hyperslices X) must be positive\n"); return 1; }

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
            Y[n2] = 0.0;
            for (int n1=0; n1<2*N1; n1+=2) { Y[n2] += (double)(X[n+n1*K]!=0.0 || X[n+n1*K+1]!=0.0); }
        }
    }

    return 0;
}


int cnt_inplace_s (float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in cnt_inplace_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in cnt_inplace_s: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in cnt_inplace_s: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in cnt_inplace_s: H (num hyperslices X) must be positive\n"); return 1; }

    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = N/(M*N1);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    float sm;

    for (int l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
    {
        for (int m=0; m<M; m++, n+=J, n2++)
        {
            sm = 0.0f;
            for (int n1=0; n1<N1; n1++) { sm += (float)(X[n+n1*K]!=0.0f); }
            X[n2] = sm;
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int cnt_inplace_d (double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in cnt_inplace_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in cnt_inplace_d: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in cnt_inplace_d: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in cnt_inplace_d: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = N/(M*N1);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    double sm;

    for (int l=0, n=0, n2=0; l<L; l++, n+=M*(N1-J))
    {
        for (int m=0; m<M; m++, n+=J, n2++)
        {
            sm = 0.0;
            for (int n1=0; n1<N1; n1++) { sm += (double)(X[n+n1*K]!=0.0); }
            X[n2] = sm;
        }
    }

    return 0;
}


int cnt_inplace_c (float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in cnt_inplace_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in cnt_inplace_c: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in cnt_inplace_c: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in cnt_inplace_c: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = N/(M*N1);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    float sm;

    for (int l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
    {
        for (int m=0; m<M; m++, n+=2*J, n2++)
        {
            sm = 0.0f;
            for (int n1=0; n1<2*N1; n1+=2) { sm += (float)(X[n+n1*K]!=0.0f || X[n+n1*K+1]!=0.0f); }
            X[n2] = sm;
        }
    }

    return 0;
}


int cnt_inplace_z (double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    if (R<1) { fprintf(stderr,"error in cnt_inplace_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in cnt_inplace_z: C (ncols X) must be positive\n"); return 1; }
    if (S<1) { fprintf(stderr,"error in cnt_inplace_z: S (num slices X) must be positive\n"); return 1; }
    if (H<1) { fprintf(stderr,"error in cnt_inplace_z: H (num hyperslices X) must be positive\n"); return 1; }

    const int RC = R*C, SH = S*H, N = RC*SH;
    const int N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int L = N/(M*N1);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    double sm;

    for (int l=0, n=0, n2=0; l<L; l++, n+=2*M*(N1-J))
    {
        for (int m=0; m<M; m++, n+=2*J, n2++)
        {
            sm = 0.0;
            for (int n1=0; n1<2*N1; n1+=2) { sm += (double)(X[n+n1*K]!=0.0 || X[n+n1*K+1]!=0.0); }
            X[n2] = sm;
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
