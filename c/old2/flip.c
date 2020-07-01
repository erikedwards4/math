//Flips X along dim (reverses order of elements).
//For d=0, this is flipud. For d=1, this is fliplr.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int flip_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim);
int flip_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim);
int flip_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim);
int flip_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim);

int flip_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim);
int flip_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim);
int flip_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim);
int flip_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim);


int flip_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim)
{
    int l, m, n1 = 0, n2;
    //struct timespec tic, toc;

    //Checks
    if (R<0) { fprintf(stderr,"error in flip_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flip_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flip_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flip_s: H (num hyperslices X) must be nonnegative\n"); return 1; }
    if (dim<0 || dim>3) { fprintf(stderr,"error in flip_s: dim must be in [0 3].\n"); return 1; }
    
    //Set ints
    //clock_gettime(CLOCK_REALTIME,&tic);
    const int L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const int M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const int N = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    // if (dim==0)
    // {
    //     if (iscolmajor) { L = C*S*H; M = R; N = 1; }
    //     else { L = 1; M = R*C*S*H; N = C*S*H; }
    // }
    // else if (dim==1)
    // {
    //     if (iscolmajor) { L = S*H; M = R*C; N = R; }
    //     else { L = R; M = C*S*H; N = S*H; }
    // }
    // else if (dim==2)
    // {
    //     if (iscolmajor) { L = H; M = R*C*S; N = R*C; }
    //     else { L = R*C; M = S*H; N = H; }
    // }
    // else if (dim==3)
    // {
    //     if (iscolmajor) { L = 1; M = R*C*S*H; N = R*C*S; }
    //     else { L = R*C*S; M = H; N = 1; }
    // }
    // else
    // {
    //     fprintf(stderr,"error in flip_inplace_s: dim must be in [0 3].\n"); return 1;
    // }

    //Flip
    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1)) { cblas_scopy(R*C*S*H,X,1,Y,1); }
    else if (N==1)
    {
        for (l=0; l<L; l++)
        {
            n2 = n1 + M - 1;
            for (m=0; m<M; m++) { Y[n2] = X[n1]; n1++; n2--; }
        }
    }
    else
    {
        for (l=0; l<L; l++)
        {
            n2 = n1 + M - N;
            for (m=0; m<M/N; m++)
            {
                cblas_scopy(N,&X[n1],1,&Y[n2],1);
                n1 += N; n2 -= N;
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int flip_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim)
{
    int l, m, n1 = 0, n2;

    //Checks
    if (R<0) { fprintf(stderr,"error in flip_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flip_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flip_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flip_d: H (num hyperslices X) must be nonnegative\n"); return 1; }
    if (dim<0 || dim>3) { fprintf(stderr,"error in flip_d: dim must be in [0 3].\n"); return 1; }

    //Set ints
    const int L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const int M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const int N = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    //Flip
    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1)) { cblas_dcopy(R*C*S*H,X,1,Y,1); }
    else if (N==1)
    {
        for (l=0; l<L; l++)
        {
            n2 = n1 + M - 1;
            for (m=0; m<M; m++) { Y[n2] = X[n1]; n1++; n2--; }
        }
    }
    else
    {
        for (l=0; l<L; l++)
        {
            n2 = n1 + M - N;
            for (m=0; m<M/N; m++)
            {
                cblas_dcopy(N,&X[n1],1,&Y[n2],1);
                n1 += N; n2 -= N;
            }
        }
    }

    return 0;
}


int flip_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim)
{
    int l, m, n1 = 0, n2;

    //Checks
    if (R<0) { fprintf(stderr,"error in flip_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flip_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flip_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flip_c: H (num hyperslices X) must be nonnegative\n"); return 1; }
    if (dim<0 || dim>3) { fprintf(stderr,"error in flip_c: dim must be in [0 3].\n"); return 1; }
    
    //Set ints
    const int L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const int M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const int N = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    //Flip
    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1)) { cblas_ccopy(R*C*S*H,X,1,Y,1); }
    else if (N==1)
    {
        for (l=0; l<L; l++)
        {
            n2 = n1 + M - 1;
            for (m=0; m<M; m++) { Y[2*n2] = X[2*n1]; Y[2*n2+1] = X[2*n1+1]; n1++; n2--; }
        }
    }
    else
    {
        for (l=0; l<L; l++)
        {
            n2 = n1 + M - N;
            for (m=0; m<M/N; m++)
            {
                cblas_ccopy(N,&X[2*n1],1,&Y[2*n2],1);
                n1 += N; n2 -= N;
            }
        }
    }

    return 0;
}


int flip_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim)
{
    int l, m, n1 = 0, n2;

    //Checks
    if (R<0) { fprintf(stderr,"error in flip_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flip_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flip_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flip_z: H (num hyperslices X) must be nonnegative\n"); return 1; }
    if (dim<0 || dim>3) { fprintf(stderr,"error in flip_z: dim must be in [0 3].\n"); return 1; }
    
    //Set ints
    const int L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const int M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const int N = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    //Flip
    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1)) { cblas_zcopy(R*C*S*H,X,1,Y,1); }
    else if (N==1)
    {
        for (l=0; l<L; l++)
        {
            n2 = n1 + M - 1;
            for (m=0; m<M; m++) { Y[2*n2] = X[2*n1]; Y[2*n2+1] = X[2*n1+1]; n1++; n2--; }
        }
    }
    else
    {
        for (l=0; l<L; l++)
        {
            n2 = n1 + M - N;
            for (m=0; m<M/N; m++)
            {
                cblas_zcopy(N,&X[2*n1],1,&Y[2*n2],1);
                n1 += N; n2 -= N;
            }
        }
    }

    return 0;
}


int flip_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim)
{
    int l, m, n1, n2;
    float x1, *X1;
    //struct timespec tic, toc;

    //Checks
    if (R<0) { fprintf(stderr,"error in flip_inplace_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flip_inplace_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flip_inplace_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flip_inplace_s: H (num hyperslices X) must be nonnegative\n"); return 1; }
    if (dim<0 || dim>3) { fprintf(stderr,"error in flip_inplace_s: dim must be in [0 3].\n"); return 1; }

    //Set ints
    //clock_gettime(CLOCK_REALTIME,&tic);
    const int L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const int M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const int N = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    //Flip
    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1)) {}
    else if (N==1)
    {
        for (l=0; l<L; l++)
        {
            n1 = l * M;
            n2 = n1 + M - 1;
            for (m=0; m<M/2; m++)
            {
                x1 = X[n1]; X[n1] = X[n2]; X[n2] = x1;
                n1++; n2--;
            }
        }
    }
    else
    {
        if (!(X1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in flip_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        for (l=0; l<L; l++)
        {
            n1 = l * M;
            n2 = n1 + M - N;
            for (m=0; m<M/N/2; m++)
            {
                //for (n=0; n<N; n++, n1++, n2--) { x1 = X[n1]; X[n1] = X[n2]; X[n2] = x1; }
                cblas_scopy(N,&X[n1],1,X1,1);
                cblas_scopy(N,&X[n2],1,&X[n1],1);
                cblas_scopy(N,X1,1,&X[n2],1);
                n1 += N; n2 -= N;
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int flip_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim)
{
    int l, m, n1, n2;
    double x1, *X1;

    //Checks
    if (R<0) { fprintf(stderr,"error in flip_inplace_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flip_inplace_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flip_inplace_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flip_inplace_d: H (num hyperslices X) must be nonnegative\n"); return 1; }
    if (dim<0 || dim>3) { fprintf(stderr,"error in flip_inplace_d: dim must be in [0 3].\n"); return 1; }

    //Set ints
    const int L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const int M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const int N = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    //Flip
    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1)) {}
    else if (N==1)
    {
        for (l=0; l<L; l++)
        {
            n1 = l * M;
            n2 = n1 + M - 1;
            for (m=0; m<M/2; m++)
            {
                x1 = X[n1]; X[n1] = X[n2]; X[n2] = x1;
                n1++; n2--;
            }
        }
    }
    else
    {
        if (!(X1=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in flip_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        for (l=0; l<L; l++)
        {
            n1 = l * M;
            n2 = n1 + M - N;
            for (m=0; m<M/N/2; m++)
            {
                cblas_dcopy(N,&X[n1],1,X1,1);
                cblas_dcopy(N,&X[n2],1,&X[n1],1);
                cblas_dcopy(N,X1,1,&X[n2],1);
                n1 += N; n2 -= N;
            }
        }
    }

    return 0;
}


int flip_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim)
{
    int l, m, n1, n2;
    float x1, *X1;

    //Checks
    if (R<0) { fprintf(stderr,"error in flip_inplace_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flip_inplace_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flip_inplace_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flip_inplace_c: H (num hyperslices X) must be nonnegative\n"); return 1; }
    if (dim<0 || dim>3) { fprintf(stderr,"error in flip_inplace_c: dim must be in [0 3].\n"); return 1; }

    //Set ints
    const int L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const int M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const int N = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    //Flip
    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1)) {}
    else if (N==1)
    {
        for (l=0; l<L; l++)
        {
            n1 = l * M;
            n2 = n1 + M - 1;
            for (m=0; m<M/2; m++)
            {
                x1 = X[2*n1]; X[2*n1] = X[2*n2]; X[2*n2] = x1;
                x1 = X[2*n1+1]; X[2*n1+1] = X[2*n2+1]; X[2*n2+1] = x1;
                n1++; n2--;
            }
        }
    }
    else
    {
        if (!(X1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in flip_inplace_c: problem with malloc. "); perror("malloc"); return 1; }
        for (l=0; l<L; l++)
        {
            n1 = l * M;
            n2 = n1 + M - N;
            for (m=0; m<M/N/2; m++)
            {
                cblas_ccopy(N,&X[2*n1],1,X1,1);
                cblas_ccopy(N,&X[2*n2],1,&X[2*n1],1);
                cblas_ccopy(N,X1,1,&X[2*n2],1);
                n1 += N; n2 -= N;
            }
        }
    }

    return 0;
}


int flip_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim)
{
    int l, m, n1, n2;
    double x1, *X1;

    //Checks
    if (R<0) { fprintf(stderr,"error in flip_inplace_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in flip_inplace_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in flip_inplace_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in flip_inplace_z: H (num hyperslices X) must be nonnegative\n"); return 1; }
    if (dim<0 || dim>3) { fprintf(stderr,"error in flip_inplace_z: dim must be in [0 3].\n"); return 1; }

    //Set ints
    const int L = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1) : ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S);
    const int M = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? R*C : (dim==2) ? R*C*S : R*C*S*H) : ((dim==0) ? R*C*S*H : (dim==1) ? C*S*H : (dim==2) ? S*H : H);
    const int N = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);

    //Flip
    if ((dim==0 && R==1) || (dim==1 && C==1) || (dim==2 && S==1) || (dim==3 && H==1)) {}
    else if (N==1)
    {
        for (l=0; l<L; l++)
        {
            n1 = l * M;
            n2 = n1 + M - 1;
            for (m=0; m<M/2; m++)
            {
                x1 = X[2*n1]; X[2*n1] = X[2*n2]; X[2*n2] = x1;
                x1 = X[2*n1+1]; X[2*n1+1] = X[2*n2+1]; X[2*n2+1] = x1;
                n1++; n2--;
            }
        }
    }
    else
    {
        if (!(X1=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in flip_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
        for (l=0; l<L; l++)
        {
            n1 = l * M;
            n2 = n1 + M - N;
            for (m=0; m<M/N/2; m++)
            {
                cblas_zcopy(N,&X[2*n1],1,X1,1);
                cblas_zcopy(N,&X[2*n2],1,&X[2*n1],1);
                cblas_zcopy(N,X1,1,&X[2*n2],1);
                n1 += N; n2 -= N;
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
