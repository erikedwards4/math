//Sorts X along dim.
//This has in-place and not-in-place versions.
//Too bad there is no LAPACKE_?laord for the in-place version.

#include <stdio.h>
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

int sort_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sort_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sort_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sort_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);

int sort_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sort_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sort_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);
int sort_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend);


int cmp_ascend_s (const void *a, const void *b)
{
	const float x1 = *(const float*)a, x2 = *(const float*)b;
	if (x1!=x1) { return 1; }
    else if (x2!=x2) { return -1; }
    else if (x1>x2) { return 1; }
    else if (x2>x1) { return -1; }
    else { return 0; }
}


int cmp_ascend_d (const void *a, const void *b)
{
	const double x1 = *(const double*)a, x2 = *(const double*)b;
	if (x1!=x1) { return 1; }
    else if (x2!=x2) { return -1; }
    else if (x1>x2) { return 1; }
    else if (x2>x1) { return -1; }
    else { return 0; }
}


int cmp_ascend_c (const void *a, const void *b)
{
	const float abs1 = sqrtf(*(const float *)a**(const float *)a + *((const float *)(a)+1)**((const float *)(a)+1));
    const float abs2 = sqrtf(*(const float *)b**(const float *)b + *((const float *)(b)+1)**((const float *)(b)+1));
	if (abs1!=abs1) { return 1; }
    else if (abs2!=abs2) { return -1; }
    else if (abs1>abs2) { return 1; }
    else if (abs2>abs1) { return -1; }
	else if (atan2f(*((const float *)(a)+1),*(const float *)a) > atan2f(*((const float *)(b)+1),*(const float *)b)) { return 1; }
    else if (atan2f(*((const float *)(b)+1),*(const float *)b) > atan2f(*((const float *)(a)+1),*(const float *)a)) { return -1; }
    else { return 0; }
}


int cmp_ascend_z (const void *a, const void *b)
{
	const double abs1 = sqrt(*(const double *)a**(const double *)a + *((const double *)(a)+1)**((const double *)(a)+1));
    const double abs2 = sqrt(*(const double *)b**(const double *)b + *((const double *)(b)+1)**((const double *)(b)+1));
	if (abs1!=abs1) { return 1; }
    else if (abs2!=abs2) { return -1; }
    else if (abs1>abs2) { return 1; }
    else if (abs2>abs1) { return -1; }
	else if (atan2(*((const double *)(a)+1),*(const double *)a) > atan2(*((const double *)(b)+1),*(const double *)b)) { return 1; }
    else if (atan2(*((const double *)(b)+1),*(const double *)b) > atan2(*((const double *)(a)+1),*(const double *)a)) { return -1; }
    else { return 0; }
}


int cmp_descend_s (const void *a, const void *b)
{
	const float x1 = *(const float*)a, x2 = *(const float*)b;
	if (x1!=x1) { return -1; }
    else if (x2!=x2) { return 1; }
    else if (x1<x2) { return 1; }
    else if (x2<x1) { return -1; }
    else { return 0; }
}


int cmp_descend_d (const void *a, const void *b)
{
	const double x1 = *(const double*)a, x2 = *(const double*)b;
	if (x1!=x1) { return -1; }
    else if (x2!=x2) { return 1; }
    else if (x1<x2) { return 1; }
    else if (x2<x1) { return -1; }
    else { return 0; }
}


int cmp_descend_c (const void *a, const void *b)
{
	const float abs1 = sqrtf(*(const float *)a**(const float *)a + *((const float *)(a)+1)**((const float *)(a)+1));
    const float abs2 = sqrtf(*(const float *)b**(const float *)b + *((const float *)(b)+1)**((const float *)(b)+1));
	if (abs1!=abs1) { return -1; }
    else if (abs2!=abs2) { return 1; }
    else if (abs1<abs2) { return 1; }
    else if (abs2<abs1) { return -1; }
	else if (atan2f(*((const float *)(a)+1),*(const float *)a) < atan2f(*((const float *)(b)+1),*(const float *)b)) { return 1; }
    else if (atan2f(*((const float *)(b)+1),*(const float *)b) < atan2f(*((const float *)(a)+1),*(const float *)a)) { return -1; }
    else { return 0; }
}


int cmp_descend_z (const void *a, const void *b)
{
	const double abs1 = sqrt(*(const double *)a**(const double *)a + *((const double *)(a)+1)**((const double *)(a)+1));
    const double abs2 = sqrt(*(const double *)b**(const double *)b + *((const double *)(b)+1)**((const double *)(b)+1));
	if (abs1!=abs1) { return -1; }
    else if (abs2!=abs2) { return 1; }
    else if (abs1<abs2) { return 1; }
    else if (abs2<abs1) { return -1; }
	else if (atan2(*((const double *)(a)+1),*(const double *)a) < atan2(*((const double *)(b)+1),*(const double *)b)) { return 1; }
    else if (atan2(*((const double *)(b)+1),*(const double *)b) < atan2(*((const double *)(a)+1),*(const double *)a)) { return -1; }
    else { return 0; }
}


int sort_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sort_s: dim must be in [0 3]\n"); return 1; }

    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t L = R*C*S*H/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const char id = (ascend) ? 'I' : 'D';
    //int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_s : cmp_descend_s;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);
    
    if (N1==1) { cblas_scopy((int)(R*C*S*H),X,1,Y,1); }
    else if (M==1 && L==1)
    {
        cblas_scopy((int)(R*C*S*H),X,1,Y,1);
        if (LAPACKE_slasrt_work(id,(int)N1,Y)) { fprintf(stderr,"error in sort_s: problem with LAPACKE function\n"); }
        //qsort(X,N1,sizeof(float),comp);
    }
    else if (K==1)
    {
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n=l*M*N1; m<M; m++, n+=J)
            {
                cblas_scopy((int)N1,&X[n],1,&Y[n],1);
                if (LAPACKE_slasrt_work(id,(int)N1,&Y[n])) { fprintf(stderr,"error in sort_s: problem with LAPACKE function\n"); }
                //qsort(&Y[n],N1,sizeof(float),comp);
            }
        }
    }
    else
    {
        float *X1;
        if (!(X1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in sort_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n=l*M*N1; m<M; m++, n+=J)
            {
                cblas_scopy((int)N1,&X[n],(int)K,X1,1);
                if (LAPACKE_slasrt_work(id,(int)N1,X1)) { fprintf(stderr,"error in sort_s: problem with LAPACKE function\n"); }
                //qsort(X1,N1,sizeof(float),comp);
                cblas_scopy((int)N1,X1,1,&Y[n],(int)K);
            }
        }
        free(X1);
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int sort_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sort_d: dim must be in [0 3]\n"); return 1; }

    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t L = R*C*S*H/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const char id = (ascend) ? 'I' : 'D';
    //int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_d : cmp_descend_d;
    
    if (N1==1) { cblas_dcopy((int)(R*C*S*H),X,1,Y,1); }
    else if (M==1 && L==1)
    {
        cblas_dcopy((int)(R*C*S*H),X,1,Y,1);
        if (LAPACKE_dlasrt(id,(int)N1,Y)) { fprintf(stderr,"error in sort_d: problem with LAPACKE function\n"); }
        //qsort(X,N1,sizeof(double),comp);
    }
    else if (K==1)
    {
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n=l*M*N1; m<M; m++, n+=J)
            {
                cblas_dcopy((int)N1,&X[n],1,&Y[n],1);
                if (LAPACKE_dlasrt(id,(int)N1,&Y[n])) { fprintf(stderr,"error in sort_d: problem with LAPACKE function\n"); }
                //qsort(&Y[n],N1,sizeof(double),comp);
            }
        }
    }
    else
    {
        double *X1;
        if (!(X1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in sort_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n=l*M*N1; m<M; m++, n+=J)
            {
                cblas_dcopy((int)N1,&X[n],(int)K,X1,1);
                if (LAPACKE_dlasrt(id,(int)N1,X1)) { fprintf(stderr,"error in sort_d: problem with LAPACKE function\n"); }
                //qsort(X1,N1,sizeof(double),comp);
                cblas_dcopy((int)N1,X1,1,&Y[n],(int)K);
            }
        }
        free(X1);
    }
    
    return 0;
}


int sort_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sort_c: dim must be in [0 3]\n"); return 1; }

    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_ccopy((int)(R*C*S*H),X,1,Y,1); }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t L = R*C*S*H/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_c : cmp_descend_c;

        float *X1;
        if (!(X1=(float *)malloc(N1*2*sizeof(float)))) { fprintf(stderr,"error in sort_c: problem with malloc. "); perror("malloc"); return 1; }
        
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n=2*l*M*N1; m<M; m++, n+=2*J)
            {
                cblas_ccopy((int)N1,&X[n],(int)K,X1,1);
                qsort(X1,N1,2*sizeof(float),comp);
                cblas_ccopy((int)N1,X1,1,&Y[n],(int)K);
            }
        }
        free(X1);
    }
    
    return 0;
}


int sort_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sort_z: dim must be in [0 3]\n"); return 1; }

    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    if (N1==1) { cblas_zcopy((int)(R*C*S*H),X,1,Y,1); }
    else
    {
        const size_t M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t L = R*C*S*H/(M*N1);
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
        int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_z : cmp_descend_z;

        double *X1;
        if (!(X1=(double *)malloc(N1*2*sizeof(double)))) { fprintf(stderr,"error in sort_z: problem with malloc. "); perror("malloc"); return 1; }
        
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n=2*l*M*N1; m<M; m++, n+=2*J)
            {
                cblas_zcopy((int)N1,&X[n],(int)K,X1,1);
                qsort(X1,N1,2*sizeof(double),comp);
                cblas_zcopy((int)N1,X1,1,&Y[n],(int)K);
            }
        }
        free(X1);
    }
    
    return 0;
}


int sort_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sort_inplace_s: dim must be in [0 3]\n"); return 1; }

    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t L = R*C*S*H/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const char id = (ascend) ? 'I' : 'D';
    //int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_s : cmp_descend_s;

    if (N1==1) {}
    else if (M==1 && L==1)
    {
        if (LAPACKE_slasrt_work(id,(int)N1,X)) { fprintf(stderr,"error in sort_inplace_s: problem with LAPACKE function\n"); }
        //qsort(X,N1,sizeof(float),comp);
    }
    else if (K==1)
    {
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n=l*M*N1; m<M; m++, n+=J)
            {
                if (LAPACKE_slasrt_work(id,(int)N1,&X[n])) { fprintf(stderr,"error in sort_inplace_s: problem with LAPACKE function\n"); }
                //qsort(&X[n],N1,sizeof(float),comp);
            }
        }
    }
    else
    {
        float *X1;
        if (!(X1=(float *)malloc(N1*sizeof(float)))) { fprintf(stderr,"error in sort_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n=l*M*N1; m<M; m++, n+=J)
            {
                cblas_scopy((int)N1,&X[n],(int)K,X1,1);
                if (LAPACKE_slasrt_work(id,(int)N1,X1)) { fprintf(stderr,"error in sort_inplace_s: problem with LAPACKE function\n"); }
                //qsort(X1,N1,sizeof(float),comp);
                cblas_scopy((int)N1,X1,1,&X[n],(int)K);
            }
        }
        free(X1);
    }

    return 0;
}


int sort_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sort_inplace_d: dim must be in [0 3]\n"); return 1; }

    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t L = R*C*S*H/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    const char id = (ascend) ? 'I' : 'D';

    if (N1==1) {}
    else if (M==1 && L==1)
    {
        if (LAPACKE_dlasrt(id,(int)N1,X)) { fprintf(stderr,"error in sort_inplace_d: problem with LAPACKE function\n"); }
    }
    else if (K==1)
    {
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n=l*M*N1; m<M; m++, n+=J)
            {
                if (LAPACKE_dlasrt(id,(int)N1,&X[n])) { fprintf(stderr,"error in sort_inplace_d: problem with LAPACKE function\n"); }
            }
        }
    }
    else
    {
        double *X1;
        if (!(X1=(double *)malloc(N1*sizeof(double)))) { fprintf(stderr,"error in sort_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n=l*M*N1; m<M; m++, n+=J)
            {
                cblas_dcopy((int)N1,&X[n],(int)K,X1,1);
                if (LAPACKE_dlasrt(id,(int)N1,X1)) { fprintf(stderr,"error in sort_inplace_d: problem with LAPACKE function\n"); }
                cblas_dcopy((int)N1,X1,1,&X[n],(int)K);
            }
        }
        free(X1);
    }

    return 0;
}


int sort_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sort_inplace_c: dim must be in [0 3]\n"); return 1; }

    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t L = R*C*S*H/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_c : cmp_descend_c;

    float *X1;
    if (!(X1=(float *)malloc(N1*2*sizeof(float)))) { fprintf(stderr,"error in sort_inplace_c: problem with malloc. "); perror("malloc"); return 1; }

    if (N1==1) {}
    else
    {
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n=2*l*M*N1; m<M; m++, n+=2*J)
            {
                cblas_ccopy((int)N1,&X[n],(int)K,X1,1);
                qsort(X1,N1,2*sizeof(float),comp);
                cblas_ccopy((int)N1,X1,1,&X[n],(int)K);
            }
        }
    }
    
    free(X1);
    return 0;
}


int sort_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const char ascend)
{
    if (dim>3) { fprintf(stderr,"error in sort_inplace_z: dim must be in [0 3]\n"); return 1; }

    const size_t N1 = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t M = (iscolmajor) ? ((dim==0) ? C*S*H : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t L = R*C*S*H/(M*N1);
    const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
    const size_t J = (iscolmajor) ? ((dim==0) ? R : (dim==1) ? 1 : (dim==2) ? 1 : 1) : ((dim==0) ? 1 : (dim==1) ? 1 : (dim==2) ? 1 : H);
    //const char id = (ascend) ? 'I' : 'D';
    int (*comp)(const void *, const void *) = (ascend) ? cmp_ascend_z : cmp_descend_z;
    
    double *X1;
    if (!(X1=(double *)malloc(N1*2*sizeof(double)))) { fprintf(stderr,"error in sort_inplace_z: problem with malloc. "); perror("malloc"); return 1; }

    if (N1==1) {}
    else
    {
        for (size_t l=0; l<L; l++)
        {
            for (size_t m=0, n=2*l*M*N1; m<M; m++, n+=2*J)
            {
                cblas_zcopy((int)N1,&X[n],(int)K,X1,1);
                qsort(X1,N1,2*sizeof(double),comp);
                cblas_zcopy((int)N1,X1,1,&X[n],(int)K);
            }
        }
    }
    
    free(X1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
