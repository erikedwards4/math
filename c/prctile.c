//Gets prctile of each row or col of X according to dim.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int cmp_ascend_s (const void *a, const void *b);
int cmp_ascend_d (const void *a, const void *b);

int prctile_s (float *Y, const float *X, const int R, const int C, const int dim, const char iscolmajor, const float p);
int prctile_d (double *Y, const double *X, const int R, const int C, const int dim, const char iscolmajor, const double p);


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


int prctile_s (float *Y, const float *X, const int R, const int C, const int dim, const char iscolmajor, const float p)
{
    int r, c, n1, n2;
    float w1, w2;
    float *X1; //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in prctile_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in prctile_s: C (ncols X) must be positive\n"); return 1; }

    //Prep interpolation weights
    w2 = (dim==0) ? (p/100.0f)*(R-1) : (p/100.0f)*(C-1);
    n1 = (int)floorf(w2); n2 = n1 + 1;
    w2 -= floorf(w2); w1 = 1.0f - w2;

    if (dim==0)
    {
        if (!(X1=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in prctile_s: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c*R],1,X1,1);
                qsort(X1,(size_t)(R),sizeof(float),cmp_ascend_s);
                Y[c] = w1*X1[n1] + w2*X1[n2];
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c],C,X1,1);
                qsort(X1,(size_t)(R),sizeof(float),cmp_ascend_s);
                Y[c] = w1*X1[n1] + w2*X1[n2];
            }
        }
    }
    else if (dim==1)
    {
        if (!(X1=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in prctile_s: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r],R,X1,1);
                qsort(X1,(size_t)(C),sizeof(float),cmp_ascend_s);
                Y[r] = w1*X1[n1] + w2*X1[n2];
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r*C],1,X1,1);
                qsort(X1,(size_t)(C),sizeof(float),cmp_ascend_s);
                Y[r] = w1*X1[n1] + w2*X1[n2];
            }
        }
    }
    else
    {
        fprintf(stderr,"error in prctile_s: dim must be 0 or 1.\n"); return 1;
    }

    free(X1);
    return 0;
}


int prctile_d (double *Y, const double *X, const int R, const int C, const int dim, const char iscolmajor, const double p)
{
    int r, c, n1, n2;
    double w1, w2;
    double *X1; //1 row or col of X

    //Checks
    if (R<1) { fprintf(stderr,"error in prctile_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in prctile_d: C (ncols X) must be positive\n"); return 1; }

    //Prep interpolation weights
    w2 = (dim==0) ? (p/100.0)*(R-1) : (p/100.0)*(C-1);
    n1 = (int)floor(w2); n2 = n1 + 1;
    w2 -= floor(w2); w1 = 1.0 - w2;

    if (dim==0)
    {
        if (!(X1=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in prctile_d: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c*R],1,X1,1);
                qsort(X1,(size_t)(R),sizeof(double),cmp_ascend_s);
                Y[c] = w1*X1[n1] + w2*X1[n2];
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c],C,X1,1);
                qsort(X1,(size_t)(R),sizeof(double),cmp_ascend_s);
                Y[c] = w1*X1[n1] + w2*X1[n2];
            }
        }
    }
    else if (dim==1)
    {
        if (!(X1=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in prctile_d: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r],R,X1,1);
                qsort(X1,(size_t)(C),sizeof(double),cmp_ascend_s);
                Y[r] = w1*X1[n1] + w2*X1[n2];
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r*C],1,X1,1);
                qsort(X1,(size_t)(C),sizeof(double),cmp_ascend_s);
                Y[r] = w1*X1[n1] + w2*X1[n2];
            }
        }
    }
    else
    {
        fprintf(stderr,"error in prctile_d: dim must be 0 or 1.\n"); return 1;
    }

    free(X1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
