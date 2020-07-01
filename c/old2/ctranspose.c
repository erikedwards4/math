//Hermitian (complex conjugate) transpose of matrix X.
//This has in-place and not-in-place versions.

//For the inplace R!=C case, there is a faster way to do this (with much more code). Perhaps consider in future.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ctranspose_s (float *Y, const float *X, const int R, const int C, const char iscolmajor);
int ctranspose_d (double *Y, const double *X, const int R, const int C, const char iscolmajor);
int ctranspose_c (float *Y, const float *X, const int R, const int C, const char iscolmajor);
int ctranspose_z (double *Y, const double *X, const int R, const int C, const char iscolmajor);

int ctranspose_inplace_s (float *X, const int R, const int C, const char iscolmajor);
int ctranspose_inplace_d (double *X, const int R, const int C, const char iscolmajor);
int ctranspose_inplace_c (float *X, const int R, const int C, const char iscolmajor);
int ctranspose_inplace_z (double *X, const int R, const int C, const char iscolmajor);


int ctranspose_s (float *Y, const float *X, const int R, const int C, const char iscolmajor)
{
    int r, c, n = 0;
    struct timespec tic, toc;
    clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (R<1) { fprintf(stderr,"error in ctranspose_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ctranspose_s: C (ncols X) must be positive\n"); return 1; }

    if (R==1 || C==1) { cblas_scopy(R*C,X,1,Y,1); }
    else if (iscolmajor)
    {
        for (c=0; c<C; c++, n+=R) { cblas_scopy(R,&X[n],1,&Y[c],C); }
    }
    else
    {
        for (r=0; r<R; r++, n+=C) { cblas_scopy(C,&X[n],1,&Y[r],R); }
    }
    clock_gettime(CLOCK_REALTIME,&toc);
    fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int ctranspose_d (double *Y, const double *X, const int R, const int C, const char iscolmajor)
{
    int r, c, n = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in ctranspose_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ctranspose_d: C (ncols X) must be positive\n"); return 1; }

    if (R==1 || C==1) { cblas_dcopy(R*C,X,1,Y,1); }
    else if (iscolmajor)
    {
        for (c=0; c<C; c++, n+=R) { cblas_dcopy(R,&X[n],1,&Y[c],C); }
    }
    else
    {
        for (r=0; r<R; r++, n+=C) { cblas_dcopy(C,&X[n],1,&Y[r],R); }
    }

    return 0;
}


int ctranspose_c (float *Y, const float *X, const int R, const int C, const char iscolmajor)
{
    int r, c, n = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in ctranspose_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ctranspose_c: C (ncols X) must be positive\n"); return 1; }

    if (R==1 || C==1) { cblas_ccopy(R*C,X,1,Y,1); }
    else if (iscolmajor)
    {
        for (c=0; c<C; c++, n+=2*R) { cblas_ccopy(R,&X[n],1,&Y[2*c],C); }
    }
    else
    {
        for (r=0; r<R; r++, n+=2*C) { cblas_ccopy(C,&X[n],1,&Y[2*r],R); }
    }
    cblas_sscal(R*C,-1.0f,&Y[1],2);

    return 0;
}


int ctranspose_z (double *Y, const double *X, const int R, const int C, const char iscolmajor)
{
    int r, c, n = 0;

    //Checks
    if (R<1) { fprintf(stderr,"error in ctranspose_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ctranspose_z: C (ncols X) must be positive\n"); return 1; }

    if (R==1 || C==1) { cblas_zcopy(R*C,X,1,Y,1); }
    else if (iscolmajor)
    {
        for (c=0; c<C; c++, n+=2*R) { cblas_zcopy(R,&X[n],1,&Y[2*c],C); }
    }
    else
    {
        for (r=0; r<R; r++, n+=2*C) { cblas_zcopy(C,&X[n],1,&Y[2*r],R); }
    }
    cblas_dscal(R*C,-1.0,&Y[1],2);

    return 0;
}


int ctranspose_inplace_s (float *X, const int R, const int C, const char iscolmajor)
{
    int r, c, n = 0;
    float *Xt;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    //Checks
    if (R<1) { fprintf(stderr,"error in ctranspose_inplace_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ctranspose_inplace_s: C (ncols X) must be positive\n"); return 1; }

    if (R==1 || C==1) {}
    else if (R==C)
    {
        for (c=0, n=1; c<C-1; c++, n+=R+1)
        {
            cblas_sswap(R-c-1,&X[n],1,&X[(c+1)*R+c],C);
        }
    }
    else
    {
        if (!(Xt=(float *)malloc((size_t)(R*C)*sizeof(float)))) { fprintf(stderr,"error in ctranspose_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy(R*C,X,1,Xt,1);
        if (iscolmajor)
        {
            for (c=0; c<C; c++, n+=R) { cblas_scopy(R,&Xt[n],1,&X[c],C); }
        }
        else
        {
            for (r=0; r<R; r++, n+=C) { cblas_scopy(C,&Xt[n],1,&X[r],R); }
        }
        free(Xt);
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int ctranspose_inplace_d (double *X, const int R, const int C, const char iscolmajor)
{
    int r, c, n = 0;
    double *Xt;

    //Checks
    if (R<1) { fprintf(stderr,"error in ctranspose_inplace_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ctranspose_inplace_d: C (ncols X) must be positive\n"); return 1; }

    if (R==1 || C==1) {}
    else if (R==C)
    {
        for (c=0, n=1; c<C-1; c++, n+=R+1)
        {
            cblas_dswap(R-c-1,&X[n],1,&X[(c+1)*R+c],C);
        }
    }
    else
    {
        if (!(Xt=(double *)malloc((size_t)(R*C)*sizeof(double)))) { fprintf(stderr,"error in ctranspose_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy(R*C,X,1,Xt,1);
        if (iscolmajor)
        {
            for (c=0; c<C; c++, n+=R) { cblas_dcopy(R,&Xt[n],1,&X[c],C); }
        }
        else
        {
            for (r=0; r<R; r++, n+=C) { cblas_dcopy(C,&Xt[n],1,&X[r],R); }
        }
        free(Xt);
    }

    return 0;
}


int ctranspose_inplace_c (float *X, const int R, const int C, const char iscolmajor)
{
    int r, c, n = 0;
    float *Xt;

    //Checks
    if (R<1) { fprintf(stderr,"error in ctranspose_inplace_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ctranspose_inplace_c: C (ncols X) must be positive\n"); return 1; }

    if (R==1 || C==1) {}
    else if (R==C)
    {
        for (c=0, n=2; c<C-1; c++, n+=2*R+2)
        {
            cblas_cswap(R-c-1,&X[n],1,&X[2*((c+1)*R+c)],C);
        }
    }
    else
    {
        if (!(Xt=(float *)malloc((size_t)(2*R*C)*sizeof(float)))) { fprintf(stderr,"error in ctranspose_inplace_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy(R*C,X,1,Xt,1);
        cblas_sscal(R*C,-1.0f,&Xt[1],2);
        if (iscolmajor)
        {
            for (c=0; c<C; c++, n+=2*R) { cblas_ccopy(R,&Xt[n],1,&X[2*c],C); }
        }
        else
        {
            for (r=0; r<R; r++, n+=2*C) { cblas_ccopy(C,&Xt[n],1,&X[2*r],R); }
        }
        free(Xt);
    }

    return 0;
}


int ctranspose_inplace_z (double *X, const int R, const int C, const char iscolmajor)
{
    int r, c, n = 0;
    double *Xt;

    //Checks
    if (R<1) { fprintf(stderr,"error in ctranspose_inplace_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in ctranspose_inplace_z: C (ncols X) must be positive\n"); return 1; }

    if (R==1 || C==1) {}
    else if (R==C)
    {
        for (c=0, n=2; c<C-1; c++, n+=2*R+2)
        {
            cblas_zswap(R-c-1,&X[n],1,&X[2*((c+1)*R+c)],C);
        }
    }
    else
    {
        if (!(Xt=(double *)malloc((size_t)(2*R*C)*sizeof(double)))) { fprintf(stderr,"error in ctranspose_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy(R*C,X,1,Xt,1);
        cblas_dscal(R*C,-1.0,&Xt[1],2);
        if (iscolmajor)
        {
            for (c=0; c<C; c++, n+=2*R) { cblas_zcopy(R,&Xt[n],1,&X[2*c],C); }
        }
        else
        {
            for (r=0; r<R; r++, n+=2*C) { cblas_zcopy(C,&Xt[n],1,&X[2*r],R); }
        }
        free(Xt);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
