//Gets ranges along dim of X.
//The range is the 100th minus the 0th percentile (max - min).

//The in-place versions still return the same Y, but modify X during processing.
//However, it turns out to be almost the identical speed for matrices.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int range_s (float *Y, const float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int range_d (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);

int range_inplace_s (float *Y, float *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);
int range_inplace_d (double *Y, double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor);


int range_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    float *X1;
    if (!(X1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in range_s: problem with malloc. "); perror("malloc"); return 1; }
    float mn, mx;

    if (N1==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        mn = mx = X[0];
        for (size_t n=1; n<N; n++)
        {
            if (X[n]<mn) { mn = X[n]; }
            else if (X[n]>mx) { mx = X[n]; }
        }
        *Y = mx - mn;
    }
    else
    {
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                cblas_scopy((int)N1,X,(int)K,X1,1);
                mn = mx = X1[0];
                for (size_t n1=1; n1<N1; n1++)
                {
                    if (X1[n1]<mn) { mn = X1[n1]; }
                    else if (X1[n1]>mx) { mx = X1[n1]; }
                }
                *Y++ = mx - mn;
            }
        }
    }

    free(X1);
    return 0;
}


int range_d (double *Y, const double *X, const size_t R, const size_t C,const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    double *X1;
    if (!(X1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in range_d: problem with malloc. "); perror("malloc"); return 1; }
    double mn, mx;

    if (N1==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        mn = mx = X[0];
        for (size_t n=1; n<N; n++)
        {
            if (X[n]<mn) { mn = X[n]; }
            else if (X[n]>mx) { mx = X[n]; }
        }
        *Y = mx - mn;
    }
    else
    {
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                cblas_dcopy((int)N1,X,(int)K,X1,1);
                mn = mx = X1[0];
                for (size_t n1=1; n1<N1; n1++)
                {
                    if (X1[n1]<mn) { mn = X1[n1]; }
                    else if (X1[n1]>mx) { mx = X1[n1]; }
                }
                *Y++ = mx - mn;
            }
        }
    }

    free(X1);
    return 0;
}


int range_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    float mn, mx;

    if (N1==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        mn = mx = X[0];
        for (size_t n=1; n<N; n++)
        {
            if (X[n]<mn) { mn = X[n]; }
            else if (X[n]>mx) { mx = X[n]; }
        }
        *Y = mx - mn;
    }
    else if (K==1)
    {
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X-=N1-J)
            {
                mn = mx = *X++;
                for (size_t n1=1; n1<N1; n1++, X++)
                {
                    if (*X<mn) { mn = *X; }
                    else if (*X>mx) { mx = *X; }
                }
                *Y++ = mx - mn;
            }
        }
    }
    else
    {
        float *X1;
        if (!(X1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in range_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                cblas_scopy((int)N1,X,(int)K,X1,1);
                mn = mx = X1[0];
                for (size_t n1=1; n1<N1; n1++)
                {
                    if (X1[n1]<mn) { mn = X1[n1]; }
                    else if (X1[n1]>mx) { mx = X1[n1]; }
                }
                *Y++ = mx - mn;
            }
        }
        free(X1);
    }

    return 0;
}


int range_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int dim, const char iscolmajor)
{
    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    double mn, mx;

    if (N1==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        mn = mx = X[0];
        for (size_t n=1; n<N; n++)
        {
            if (X[n]<mn) { mn = X[n]; }
            else if (X[n]>mx) { mx = X[n]; }
        }
        *Y = mx - mn;
    }
    else if (K==1)
    {
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X-=N1-J)
            {
                mn = mx = *X++;
                for (size_t n1=1; n1<N1; n1++, X++)
                {
                    if (*X<mn) { mn = *X; }
                    else if (*X>mx) { mx = *X; }
                }
                *Y++ = mx - mn;
            }
        }
    }
    else
    {
        double *X1;
        if (!(X1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in range_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                cblas_dcopy((int)N1,X,(int)K,X1,1);
                mn = mx = X1[0];
                for (size_t n1=1; n1<N1; n1++)
                {
                    if (X1[n1]<mn) { mn = X1[n1]; }
                    else if (X1[n1]>mx) { mx = X1[n1]; }
                }
                *Y++ = mx - mn;
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
