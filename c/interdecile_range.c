//Gets interdecile ranges along dim of X.
//The IDR is the 90th minus the 10th percentile.

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

int interdecile_range_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int interdecile_range_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int interdecile_range_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int interdecile_range_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int interdecile_range_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in interdecile_range_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    //Prep interpolation
    const float p1 = 0.1f*(N1-1), p2 = 0.9f*(N1-1);
    const size_t n1 = (size_t)floorf(p1), n2 = n1 + 1;
    const size_t n3 = (size_t)floorf(p2), n4 = n3 + 1;
    const float w2 = p1 - floorf(p1), w1 = 1.0f - w2;
    const float w4 = p2 - floorf(p2), w3 = 1.0f - w4;

    float *X1;
    if (!(X1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in interdecile_range_s: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N1==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        cblas_scopy((int)N1,X,1,X1,1);
        if (LAPACKE_slasrt_work('I',(int)N,X1)) { fprintf(stderr,"error in interdecile_range_s: problem with LAPACKE function\n"); }
        *Y = w3*X1[n3] + w4*X1[n4] - w1*X1[n1] - w2*X1[n2];
    }
    else
    {
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                cblas_scopy((int)N1,X,(int)K,X1,1);
                if (LAPACKE_slasrt_work('I',(int)N1,X1)) { fprintf(stderr,"error in interdecile_range_s: problem with LAPACKE function\n"); }
                *Y++ = w3*X1[n3] + w4*X1[n4] - w1*X1[n1] - w2*X1[n2];
            }
        }
    }

    free(X1);
    return 0;
}


int interdecile_range_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in interdecile_range_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    //Prep interpolation
    const double p1 = 0.1*(N1-1), p2 = 0.9*(N1-1);
    const size_t n1 = (size_t)floor(p1), n2 = n1 + 1;
    const size_t n3 = (size_t)floor(p2), n4 = n3 + 1;
    const double w2 = p1 - floor(p1), w1 = 1.0 - w2;
    const double w4 = p2 - floor(p2), w3 = 1.0 - w4;

    double *X1;
    if (!(X1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in interdecile_range_d: problem with malloc. "); perror("malloc"); return 1; }
    
    if (N1==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        cblas_dcopy((int)N1,X,1,X1,1);
        if (LAPACKE_dlasrt_work('I',(int)N,X1)) { fprintf(stderr,"error in interdecile_range_d: problem with LAPACKE function\n"); }
        *Y = w3*X1[n3] + w4*X1[n4] - w1*X1[n1] - w2*X1[n2];
    }
    else
    {
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                cblas_dcopy((int)N1,X,(int)K,X1,1);
                if (LAPACKE_dlasrt_work('I',(int)N1,X1)) { fprintf(stderr,"error in interdecile_range_d: problem with LAPACKE function\n"); }
                *Y++ = w3*X1[n3] + w4*X1[n4] - w1*X1[n1] - w2*X1[n2];
            }
        }
    }

    free(X1);
    return 0;
}


int interdecile_range_inplace_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in interdecile_range_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    //Prep interpolation
    const float p1 = 0.1f*(N1-1), p2 = 0.9f*(N1-1);
    const size_t n1 = (size_t)floorf(p1), n2 = n1 + 1;
    const size_t n3 = (size_t)floorf(p2), n4 = n3 + 1;
    const float w2 = p1 - floorf(p1), w1 = 1.0f - w2;
    const float w4 = p2 - floorf(p2), w3 = 1.0f - w4;

    if (N1==1) { cblas_scopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        if (LAPACKE_slasrt_work('I',(int)N,X)) { fprintf(stderr,"error in interdecile_range_inplace_s: problem with LAPACKE function\n"); }
        *Y = w3*X[n3] + w4*X[n4] - w1*X[n1] - w2*X[n2];
    }
    else if (K==1)
    {
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J-n3)
            {
                if (LAPACKE_slasrt_work('I',(int)N1,X)) { fprintf(stderr,"error in interdecile_range_inplace_s: problem with LAPACKE function\n"); }
                X += n1; *Y = -w1**X - w2**(X+1);
                X += n3 - n1; *Y++ += w3**X + w4**(X+1);
            }
        }
    }
    else
    {
        float *X1;
        if (!(X1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in interdecile_range_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                cblas_scopy((int)N1,X,(int)K,X1,1);
                if (LAPACKE_slasrt_work('I',(int)N1,X1)) { fprintf(stderr,"error in interdecile_range_inplace_s: problem with LAPACKE function\n"); }
                *Y++ = w3*X1[n3] + w4*X1[n4] - w1*X1[n1] - w2*X1[n2];
            }
        }
        free(X1);
    }

    return 0;
}


int interdecile_range_inplace_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in interdecile_range_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t RC = R*C, SH = S*H, N = RC*SH;
    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*SH : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t L = N/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? RC : RC*S) : ((dim==0) ? C*SH : (dim==1) ? SH : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);

    //Prep interpolation
    const double p1 = 0.1*(N1-1), p2 = 0.9*(N1-1);
    const size_t n1 = (size_t)floor(p1), n2 = n1 + 1;
    const size_t n3 = (size_t)floor(p2), n4 = n3 + 1;
    const double w2 = p1 - floor(p1), w1 = 1.0 - w2;
    const double w4 = p2 - floor(p2), w3 = 1.0 - w4;

    if (N1==1) { cblas_dcopy((int)N,X,1,Y,1); }
    else if (N1==N)
    {
        if (LAPACKE_dlasrt_work('I',(int)N,X)) { fprintf(stderr,"error in interdecile_range_inplace_d: problem with LAPACKE function\n"); }
        *Y = w3*X[n3] + w4*X[n4] - w1*X[n1] - w2*X[n2];
    }
    else if (K==1)
    {
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J-n3)
            {
                if (LAPACKE_dlasrt_work('I',(int)N1,X)) { fprintf(stderr,"error in interdecile_range_inplace_d: problem with LAPACKE function\n"); }
                X += n1; *Y = -w1**X - w2**(X+1);
                X += n3 - n1; *Y++ += w3**X + w4**(X+1);
            }
        }
    }
    else
    {
        double *X1;
        if (!(X1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in interdecile_range_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++, X+=M*(N1-J))
        {
            for (size_t m=0; m<M; m++, X+=J)
            {
                cblas_dcopy((int)N1,X,(int)K,X1,1);
                if (LAPACKE_dlasrt_work('I',(int)N1,X1)) { fprintf(stderr,"error in interdecile_range_inplace_d: problem with LAPACKE function\n"); }
                *Y++ = w3*X1[n3] + w4*X1[n4] - w1*X1[n1] - w2*X1[n2];
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
