//Sorts X along dim.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cmp_ascend_s (const void *a, const void *b);
int cmp_ascend_d (const void *a, const void *b);
int cmp_ascend_c (const void *a, const void *b);
int cmp_ascend_z (const void *a, const void *b);

int cmp_descend_s (const void *a, const void *b);
int cmp_descend_d (const void *a, const void *b);
int cmp_descend_c (const void *a, const void *b);
int cmp_descend_z (const void *a, const void *b);

int sort_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend);
int sort_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend);
int sort_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend);
int sort_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend);

int sort_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend);
int sort_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend);
int sort_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend);
int sort_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend);


int cmp_ascend_s (const void *a, const void *b)
{
	float x1 = *(const float*)a, x2 = *(const float*)b;
	if (x1!=x1) { return 1; }
    else if (x2!=x2) { return -1; }
    else if (x1>x2) { return 1; }
    else if (x2>x1) { return -1; }
    else { return 0; }
}


int cmp_ascend_d (const void *a, const void *b)
{
	double x1 = *(const double*)a, x2 = *(const double*)b;
	if (x1!=x1) { return 1; }
    else if (x2!=x2) { return -1; }
    else if (x1>x2) { return 1; }
    else if (x2>x1) { return -1; }
    else { return 0; }
}


int cmp_ascend_c (const void *a, const void *b)
{
    float x1[2], x2[2], abs1, abs2;
    memcpy(x1,a,2*sizeof(float));
    memcpy(x2,b,2*sizeof(float));
    abs1 = sqrtf(x1[0]*x1[0] + x1[1]*x1[1]);
    abs2 = sqrtf(x2[0]*x2[0] + x2[1]*x2[1]);
	if (abs1!=abs1) { return 1; }
    else if (abs2!=abs2) { return -1; }
    else if (abs1>abs2) { return 1; }
    else if (abs2>abs1) { return -1; }
	else if (atan2f(x1[1],x1[0])>atan2f(x2[1],x2[0])) { return 1; }
    else if (atan2f(x2[1],x2[0])>atan2f(x1[1],x1[0])) { return -1; }
    else { return 0; }
}


int cmp_ascend_z (const void *a, const void *b)
{
    double x1[2], x2[2], abs1, abs2;
    memcpy(x1,a,2*sizeof(double));
    memcpy(x2,b,2*sizeof(double));
    abs1 = sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
    abs2 = sqrt(x2[0]*x2[0] + x2[1]*x2[1]);
	if (abs1!=abs1) { return 1; }
    else if (abs2!=abs2) { return -1; }
    else if (abs1>abs2) { return 1; }
    else if (abs2>abs1) { return -1; }
	else if (atan2(x1[1],x1[0])>atan2(x2[1],x2[0])) { return 1; }
    else if (atan2(x2[1],x2[0])>atan2(x1[1],x1[0])) { return -1; }
    else { return 0; }
}


int cmp_descend_s (const void *a, const void *b)
{
	float x1 = *(const float*)a, x2 = *(const float*)b;
	if (x1!=x1) { return -1; }
    else if (x2!=x2) { return 1; }
    else if (x1<x2) { return 1; }
    else if (x2<x1) { return -1; }
    else { return 0; }
}


int cmp_descend_d (const void *a, const void *b)
{
	double x1 = *(const double*)a, x2 = *(const double*)b;
	if (x1!=x1) { return -1; }
    else if (x2!=x2) { return 1; }
    else if (x1<x2) { return 1; }
    else if (x2<x1) { return -1; }
    else { return 0; }
}


int cmp_descend_c (const void *a, const void *b)
{
	float x1[2], x2[2], abs1, abs2;
    memcpy(x1,a,2*sizeof(float));
    memcpy(x2,b,2*sizeof(float));
    abs1 = sqrtf(x1[0]*x1[0] + x1[1]*x1[1]);
    abs2 = sqrtf(x2[0]*x2[0] + x2[1]*x2[1]);
	if (abs1!=abs1) { return -1; }
    else if (abs2!=abs2) { return 1; }
    else if (abs1<abs2) { return 1; }
    else if (abs2<abs1) { return -1; }
	else if (atan2f(x1[1],x1[0])<atan2f(x2[1],x2[0])) { return 1; }
    else if (atan2f(x2[1],x2[0])<atan2f(x1[1],x1[0])) { return -1; }
    else { return 0; }
}


int cmp_descend_z (const void *a, const void *b)
{
	double x1[2], x2[2], abs1, abs2;
    memcpy(x1,a,2*sizeof(double));
    memcpy(x2,b,2*sizeof(double));
    abs1 = sqrt(x1[0]*x1[0] + x1[1]*x1[1]);
    abs2 = sqrt(x2[0]*x2[0] + x2[1]*x2[1]);
	if (abs1!=abs1) { return -1; }
    else if (abs2!=abs2) { return 1; }
    else if (abs1<abs2) { return 1; }
    else if (abs2<abs1) { return -1; }
	else if (atan2(x1[1],x1[0])<atan2(x2[1],x2[0])) { return 1; }
    else if (atan2(x2[1],x2[0])<atan2(x1[1],x1[0])) { return -1; }
    else { return 0; }
}


int sort_s (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    int l, m, n;
    float *X1;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (R<0) { fprintf(stderr,"error in sort_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    //Allocate
    if (!(X1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in sort_s: problem with malloc. "); perror("malloc"); return 1; }

    //Sort
    if (N==1) { cblas_scopy(R*C*S*H,X,1,Y,1); }
    else if (ascend)
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=l*M*N; m<M; m++, n+=J)
            {
                cblas_scopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),sizeof(float),cmp_ascend_s);
                cblas_scopy(N,X1,1,&Y[n],I);
            }
        }
    }
    else
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=l*M*N; m<M; m++, n+=J)
            {
                cblas_scopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),sizeof(float),cmp_descend_s);
                cblas_scopy(N,X1,1,&Y[n],I);
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int sort_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    int l, m, n;
    double *X1;

    //Checks
    if (R<0) { fprintf(stderr,"error in sort_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    //Allocate
    if (!(X1=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in sort_d: problem with malloc. "); perror("malloc"); return 1; }

    //Sort
    if (N==1) { cblas_dcopy(R*C*S*H,X,1,Y,1); }
    else if (ascend)
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=l*M*N; m<M; m++, n+=J)
            {
                cblas_dcopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),sizeof(double),cmp_ascend_d);
                cblas_dcopy(N,X1,1,&Y[n],I);
            }
        }
    }
    else
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=l*M*N; m<M; m++, n+=J)
            {
                cblas_dcopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),sizeof(double),cmp_descend_d);
                cblas_dcopy(N,X1,1,&Y[n],I);
            }
        }
    }
    
    return 0;
}


int sort_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    int l, m, n;
    float *X1;

    //Checks
    if (R<0) { fprintf(stderr,"error in sort_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    //Allocate
    if (!(X1=(float *)malloc((size_t)N*2*sizeof(float)))) { fprintf(stderr,"error in sort_c: problem with malloc. "); perror("malloc"); return 1; }

    //Sort
    if (N==1) { cblas_ccopy(R*C*S*H,X,1,Y,1); }
    else if (ascend)
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=2*l*M*N; m<M; m++, n+=2*J)
            {
                cblas_ccopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),2*sizeof(float),cmp_ascend_c);
                cblas_ccopy(N,X1,1,&Y[n],I);
            }
        }
    }
    else
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=2*l*M*N; m<M; m++, n+=2*J)
            {
                cblas_ccopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),2*sizeof(float),cmp_descend_c);
                cblas_ccopy(N,X1,1,&Y[n],I);
            }
        }
    }
    
    return 0;
}


int sort_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    int l, m, n;
    double *X1;

    //Checks
    if (R<0) { fprintf(stderr,"error in sort_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    //Allocate
    if (!(X1=(double *)malloc((size_t)N*2*sizeof(double)))) { fprintf(stderr,"error in sort_z: problem with malloc. "); perror("malloc"); return 1; }

    //Sort
    if (N==1) { cblas_zcopy(R*C*S*H,X,1,Y,1); }
    else if (ascend)
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=2*l*M*N; m<M; m++, n+=2*J)
            {
                cblas_zcopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),2*sizeof(double),cmp_ascend_z);
                cblas_zcopy(N,X1,1,&Y[n],I);
            }
        }
    }
    else
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=2*l*M*N; m<M; m++, n+=2*J)
            {
                cblas_zcopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),2*sizeof(double),cmp_descend_z);
                cblas_zcopy(N,X1,1,&Y[n],I);
            }
        }
    }
    
    return 0;
}


int sort_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    int l, m, n;
    float *X1;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (R<0) { fprintf(stderr,"error in sort_inplace_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_inplace_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_inplace_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_inplace_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    //Allocate
    if (!(X1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in sort_s: problem with malloc. "); perror("malloc"); return 1; }

    //Sort
    if (N==1) {}
    else if (ascend)
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=l*M*N; m<M; m++, n+=J)
            {
                cblas_scopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),sizeof(float),cmp_ascend_s);
                cblas_scopy(N,X1,1,&X[n],I);
            }
        }
    }
    else
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=l*M*N; m<M; m++, n+=J)
            {
                cblas_scopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),sizeof(float),cmp_descend_s);
                cblas_scopy(N,X1,1,&X[n],I);
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int sort_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    int l, m, n;
    double *X1;

    //Checks
    if (R<0) { fprintf(stderr,"error in sort_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    //Allocate
    if (!(X1=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in sort_d: problem with malloc. "); perror("malloc"); return 1; }

    //Sort
    if (N==1) {}
    else if (ascend)
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=2*l*M*N; m<M; m++, n+=2*J)
            {
                cblas_dcopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),sizeof(double),cmp_ascend_d);
                cblas_dcopy(N,X1,1,&X[n],I);
            }
        }
    }
    else
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=2*l*M*N; m<M; m++, n+=2*J)
            {
                cblas_dcopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),sizeof(double),cmp_descend_d);
                cblas_dcopy(N,X1,1,&X[n],I);
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int sort_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    int l, m, n;
    float *X1;

    //Checks
    if (R<0) { fprintf(stderr,"error in sort_inplace_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_inplace_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_inplace_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_inplace_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    //Allocate
    if (!(X1=(float *)malloc((size_t)N*2*sizeof(float)))) { fprintf(stderr,"error in sort_inplace_c: problem with malloc. "); perror("malloc"); return 1; }

    //Sort
    if (N==1) {}
    else if (ascend)
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=2*l*M*N; m<M; m++, n+=2*J)
            {
                cblas_ccopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),2*sizeof(float),cmp_ascend_c);
                cblas_ccopy(N,X1,1,&X[n],I);
            }
        }
    }
    else
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=2*l*M*N; m<M; m++, n+=2*J)
            {
                cblas_ccopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),2*sizeof(float),cmp_descend_c);
                cblas_ccopy(N,X1,1,&X[n],I);
            }
        }
    }
    
    return 0;
}


int sort_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    int l, m, n;
    double *X1;

    //Checks
    if (R<0) { fprintf(stderr,"error in sort_inplace_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_inplace_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_inplace_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_inplace_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    //Set ints
    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int I = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    
    //Allocate
    if (!(X1=(double *)malloc((size_t)N*2*sizeof(double)))) { fprintf(stderr,"error in sort_inplace_z: problem with malloc. "); perror("malloc"); return 1; }

    //Sort
    if (N==1) {}
    else if (ascend)
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=2*l*M*N; m<M; m++, n+=2*J)
            {
                cblas_zcopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),2*sizeof(double),cmp_ascend_z);
                cblas_zcopy(N,X1,1,&X[n],I);
            }
        }
    }
    else
    {
        for (l=0; l<L; l++)
        {
            for (m=0, n=2*l*M*N; m<M; m++, n+=2*J)
            {
                cblas_zcopy(N,&X[n],I,X1,1);
                qsort(X1,(size_t)(N),2*sizeof(double),cmp_descend_z);
                cblas_zcopy(N,X1,1,&X[n],I);
            }
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
