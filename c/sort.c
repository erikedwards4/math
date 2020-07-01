//Sorts X along dim.
//This has in-place and not-in-place versions.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include <time.h>

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
    memcpy(x1,a,2*sizeof(float)); memcpy(x2,b,2*sizeof(float));
    abs1 = sqrtf(x1[0]*x1[0] + x1[1]*x1[1]); abs2 = sqrtf(x2[0]*x2[0] + x2[1]*x2[1]);
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
    memcpy(x1,a,2*sizeof(double)); memcpy(x2,b,2*sizeof(double));
    abs1 = sqrt(x1[0]*x1[0] + x1[1]*x1[1]); abs2 = sqrt(x2[0]*x2[0] + x2[1]*x2[1]);
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
    memcpy(x1,a,2*sizeof(float)); memcpy(x2,b,2*sizeof(float));
    abs1 = sqrtf(x1[0]*x1[0] + x1[1]*x1[1]); abs2 = sqrtf(x2[0]*x2[0] + x2[1]*x2[1]);
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
    memcpy(x1,a,2*sizeof(double)); memcpy(x2,b,2*sizeof(double));
    abs1 = sqrt(x1[0]*x1[0] + x1[1]*x1[1]); abs2 = sqrt(x2[0]*x2[0] + x2[1]*x2[1]);
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
    if (R<0) { fprintf(stderr,"error in sort_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const char id = (ascend) ? 'I' : 'D';
    //int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_s : cmp_descend_s;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    
    //Sort
    if (N==1) { cblas_scopy(R*C*S*H,X,1,Y,1); }
    else if (M==1 && L==1)
    {
        cblas_scopy(R*C*S*H,X,1,Y,1);
        if (LAPACKE_slasrt_work(id,N,Y)) { fprintf(stderr,"error in sort_s: problem with LAPACKE function\n"); }
        //qsort(X,(size_t)(N),sizeof(float),comp);
    }
    else if (K==1)
    {
        int n;
        for (int l=0; l<L; l++)
        {
            n = l*M*N;
            for (int m=0; m<M; m++, n+=J)
            {
                cblas_scopy(N,&X[n],1,&Y[n],1);
                if (LAPACKE_slasrt_work(id,N,&Y[n])) { fprintf(stderr,"error in sort_s: problem with LAPACKE function\n"); }
                //qsort(&Y[n],(size_t)(N),sizeof(float),comp);
            }
        }
    }
    else
    {
        int n;
        float *X1;
        if (!(X1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in sort_s: problem with malloc. "); perror("malloc"); return 1; }
        for (int l=0; l<L; l++)
        {
            n = l*M*N;
            for (int m=0; m<M; m++, n+=J)
            {
                cblas_scopy(N,&X[n],K,X1,1);
                if (LAPACKE_slasrt_work(id,N,X1)) { fprintf(stderr,"error in sort_s: problem with LAPACKE function\n"); }
                //qsort(X1,(size_t)(N),sizeof(float),comp);
                cblas_scopy(N,X1,1,&Y[n],K);
            }
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int sort_d (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    if (R<0) { fprintf(stderr,"error in sort_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const char id = (ascend) ? 'I' : 'D';
    //int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_d : cmp_descend_d;
    
    //Sort
    if (N==1) { cblas_dcopy(R*C*S*H,X,1,Y,1); }
    else if (M==1 && L==1)
    {
        cblas_dcopy(R*C*S*H,X,1,Y,1);
        if (LAPACKE_dlasrt(id,N,Y)) { fprintf(stderr,"error in sort_d: problem with LAPACKE function\n"); }
        //qsort(X,(size_t)(N),sizeof(double),comp);
    }
    else if (K==1)
    {
        int n;
        for (int l=0; l<L; l++)
        {
            n = l*M*N;
            for (int m=0; m<M; m++, n+=J)
            {
                cblas_dcopy(N,&X[n],1,&Y[n],1);
                if (LAPACKE_dlasrt(id,N,&Y[n])) { fprintf(stderr,"error in sort_d: problem with LAPACKE function\n"); }
                //qsort(&Y[n],(size_t)(N),sizeof(double),comp);
            }
        }
    }
    else
    {
        int n;
        double *X1;
        if (!(X1=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in sort_d: problem with malloc. "); perror("malloc"); return 1; }
        for (int l=0; l<L; l++)
        {
            n = l*M*N;
            for (int m=0; m<M; m++, n+=J)
            {
                cblas_dcopy(N,&X[n],K,X1,1);
                if (LAPACKE_dlasrt(id,N,X1)) { fprintf(stderr,"error in sort_d: problem with LAPACKE function\n"); }
                //qsort(X1,(size_t)(N),sizeof(double),comp);
                cblas_dcopy(N,X1,1,&Y[n],K);
            }
        }
    }
    
    return 0;
}


int sort_c (float *Y, const float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    if (R<0) { fprintf(stderr,"error in sort_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==1) { cblas_ccopy(R*C*S*H,X,1,Y,1); }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const int L = R*C*S*H/(M*N);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_c : cmp_descend_c;
        int n;

        float *X1;
        if (!(X1=(float *)malloc((size_t)N*2*sizeof(float)))) { fprintf(stderr,"error in sort_c: problem with malloc. "); perror("malloc"); return 1; }
        
        for (int l=0; l<L; l++)
        {
            n = 2*l*M*N;
            for (int m=0; m<M; m++, n+=2*J)
            {
                cblas_ccopy(N,&X[n],K,X1,1);
                qsort(X1,(size_t)(N),2*sizeof(float),comp);
                cblas_ccopy(N,X1,1,&Y[n],K);
            }
        }
    }
    
    return 0;
}


int sort_z (double *Y, const double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    if (R<0) { fprintf(stderr,"error in sort_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N==1) { cblas_zcopy(R*C*S*H,X,1,Y,1); }
    else
    {
        const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const int L = R*C*S*H/(M*N);
        const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_z : cmp_descend_z;
        int n;

        double *X1;
        if (!(X1=(double *)malloc((size_t)N*2*sizeof(double)))) { fprintf(stderr,"error in sort_z: problem with malloc. "); perror("malloc"); return 1; }
        
        for (int l=0; l<L; l++)
        {
            n = 2*l*M*N;
            for (int m=0; m<M; m++, n+=2*J)
            {
                cblas_zcopy(N,&X[n],K,X1,1);
                qsort(X1,(size_t)(N),2*sizeof(double),comp);
                cblas_zcopy(N,X1,1,&Y[n],K);
            }
        }
    }
    
    return 0;
}


int sort_inplace_s (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    if (R<0) { fprintf(stderr,"error in sort_inplace_s: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_inplace_s: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_inplace_s: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_inplace_s: H (num hyperslices X) must be nonnegative\n"); return 1; }

    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const char id = (ascend) ? 'I' : 'D';
    //int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_s : cmp_descend_s;

    if (N==1) {}
    else if (M==1 && L==1)
    {
        if (LAPACKE_slasrt_work(id,N,X)) { fprintf(stderr,"error in sort_inplace_s: problem with LAPACKE function\n"); }
        //qsort(X,(size_t)(N),sizeof(float),comp);
    }
    else if (K==1)
    {
        int n;
        for (int l=0; l<L; l++)
        {
            n = l*M*N;
            for (int m=0; m<M; m++, n+=J)
            {
                if (LAPACKE_slasrt_work(id,N,&X[n])) { fprintf(stderr,"error in sort_inplace_s: problem with LAPACKE function\n"); }
                //qsort(&X[n],(size_t)(N),sizeof(float),comp);
            }
        }
    }
    else
    {
        int n;
        float *X1;
        if (!(X1=(float *)malloc((size_t)N*sizeof(float)))) { fprintf(stderr,"error in sort_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        for (int l=0; l<L; l++)
        {
            n = l*M*N;
            for (int m=0; m<M; m++, n+=J)
            {
                cblas_scopy(N,&X[n],K,X1,1);
                if (LAPACKE_slasrt_work(id,N,X1)) { fprintf(stderr,"error in sort_inplace_s: problem with LAPACKE function\n"); }
                //qsort(X1,(size_t)(N),sizeof(float),comp);
                cblas_scopy(N,X1,1,&X[n],K);
            }
        }
    }

    return 0;
}


int sort_inplace_d (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    if (R<0) { fprintf(stderr,"error in sort_inplace_d: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_inplace_d: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_inplace_d: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_inplace_d: H (num hyperslices X) must be nonnegative\n"); return 1; }

    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const char id = (ascend) ? 'I' : 'D';

    if (N==1) {}
    else if (M==1 && L==1)
    {
        if (LAPACKE_dlasrt(id,N,X)) { fprintf(stderr,"error in sort_inplace_d: problem with LAPACKE function\n"); }
    }
    else if (K==1)
    {
        int n;
        for (int l=0; l<L; l++)
        {
            n = l*M*N;
            for (int m=0; m<M; m++, n+=J)
            {
                if (LAPACKE_dlasrt(id,N,&X[n])) { fprintf(stderr,"error in sort_inplace_d: problem with LAPACKE function\n"); }
            }
        }
    }
    else
    {
        int n;
        double *X1;
        if (!(X1=(double *)malloc((size_t)N*sizeof(double)))) { fprintf(stderr,"error in sort_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        for (int l=0; l<L; l++)
        {
            n = l*M*N;
            for (int m=0; m<M; m++, n+=J)
            {
                cblas_dcopy(N,&X[n],K,X1,1);
                if (LAPACKE_dlasrt(id,N,X1)) { fprintf(stderr,"error in sort_inplace_d: problem with LAPACKE function\n"); }
                cblas_dcopy(N,X1,1,&X[n],K);
            }
        }
    }

    return 0;
}


int sort_inplace_c (float *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    if (R<0) { fprintf(stderr,"error in sort_inplace_c: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_inplace_c: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_inplace_c: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_inplace_c: H (num hyperslices X) must be nonnegative\n"); return 1; }

    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_c : cmp_descend_c;

    float *X1;
    if (!(X1=(float *)malloc((size_t)N*2*sizeof(float)))) { fprintf(stderr,"error in sort_inplace_c: problem with malloc. "); perror("malloc"); return 1; }

    if (N==1) {}
    else
    {
        int n;
        for (int l=0; l<L; l++)
        {
            n = 2*l*M*N;
            for (int m=0; m<M; m++, n+=2*J)
            {
                cblas_ccopy(N,&X[n],K,X1,1);
                qsort(X1,(size_t)(N),2*sizeof(float),comp);
                cblas_ccopy(N,X1,1,&X[n],K);
            }
        }
    }
    
    return 0;
}


int sort_inplace_z (double *X, const int R, const int C, const int S, const int H, const char iscolmajor, const int dim, const char ascend)
{
    if (R<0) { fprintf(stderr,"error in sort_inplace_z: R (num rows X) must be nonnegative\n"); return 1; }
    if (C<0) { fprintf(stderr,"error in sort_inplace_z: C (num cols X) must be nonnegative\n"); return 1; }
    if (S<0) { fprintf(stderr,"error in sort_inplace_z: S (num slices X) must be nonnegative\n"); return 1; }
    if (H<0) { fprintf(stderr,"error in sort_inplace_z: H (num hyperslices X) must be nonnegative\n"); return 1; }

    const int N = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const int M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int L = R*C*S*H/(M*N);
    const int K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const int J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    //const char id = (ascend) ? 'I' : 'D';
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_z : cmp_descend_z;
    
    double *X1;
    if (!(X1=(double *)malloc((size_t)N*2*sizeof(double)))) { fprintf(stderr,"error in sort_inplace_z: problem with malloc. "); perror("malloc"); return 1; }

    if (N==1) {}
    else
    {
        int n;
        for (int l=0; l<L; l++)
        {
            n = 2*l*M*N;
            for (int m=0; m<M; m++, n+=2*J)
            {
                cblas_zcopy(N,&X[n],K,X1,1);
                qsort(X1,(size_t)(N),2*sizeof(double),comp);
                cblas_zcopy(N,X1,1,&X[n],K);
            }
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
