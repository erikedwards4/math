//Non-hermitian transpose of matrix X.
//This has in-place and not-in-place versions.

//For the inplace R!=C case, there may be a faster way to do this (with much more code).
//However, it is already ~same speed as Armadillo inplace_strans.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int transpose_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor);
int transpose_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor);
int transpose_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor);
int transpose_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor);

int transpose_inplace_s (float *X, const size_t R, const size_t C, const char iscolmajor);
int transpose_inplace_d (double *X, const size_t R, const size_t C, const char iscolmajor);
int transpose_inplace_c (float *X, const size_t R, const size_t C, const char iscolmajor);
int transpose_inplace_z (double *X, const size_t R, const size_t C, const char iscolmajor);


int transpose_s (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor)
{
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if (R==1 || C==1) { cblas_scopy((int)(R*C),X,1,Y,1); }
    else if (iscolmajor)
    {
        for (size_t c=0; c<C; ++c, X+=R, ++Y) { cblas_scopy((int)(R),X,1,Y,(int)C); }
    }
    else
    {
        for (size_t r=0; r<R; ++r, X+=C, ++Y) { cblas_scopy((int)(C),X,1,Y,(int)R); }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int transpose_d (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R==1 || C==1) { cblas_dcopy((int)(R*C),X,1,Y,1); }
    else if (iscolmajor)
    {
        for (size_t c=0; c<C; ++c, X+=R, ++Y) { cblas_dcopy((int)(R),X,1,Y,(int)C); }
    }
    else
    {
        for (size_t r=0; r<R; ++r, X+=C, ++Y) { cblas_dcopy((int)(C),X,1,Y,(int)R); }
    }

    return 0;
}


int transpose_c (float *Y, const float *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R==1 || C==1) { cblas_ccopy((int)(R*C),X,1,Y,1); }
    else if (iscolmajor)
    {
        for (size_t c=0; c<C; c+=2, X+=2*R, Y+=2) { cblas_ccopy((int)(R),X,1,Y,(int)C); }
    }
    else
    {
        for (size_t r=0; r<R; r+=2, X+=2*C, Y+=2) { cblas_ccopy((int)(C),X,1,Y,(int)R); }
    }

    return 0;
}


int transpose_z (double *Y, const double *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R==1 || C==1) { cblas_zcopy((int)(R*C),X,1,Y,1); }
    else if (iscolmajor)
    {
        for (size_t c=0; c<C; c+=2, X+=2*R, Y+=2) { cblas_zcopy((int)(R),X,1,Y,(int)C); }
    }
    else
    {
        for (size_t r=0; r<R; r+=2, X+=2*C, Y+=2) { cblas_zcopy((int)(C),X,1,Y,(int)R); }
    }

    return 0;
}


int transpose_inplace_s (float *X, const size_t R, const size_t C, const char iscolmajor)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    if (R==1 || C==1) {}
    else if (R==C)
    {
        for (size_t c=0, n=1; c<C-1; ++c, n+=R+1)
        {
            cblas_sswap((int)(R-c-1),&X[n],1,&X[(c+1)*R+c],(int)C);
        }
    }
    else
    {
        float *Xt;
        if (!(Xt=(float *)malloc(R*C*sizeof(float)))) { fprintf(stderr,"error in transpose_inplace_s: problem with malloc. "); perror("malloc"); return 1; }
        cblas_scopy((int)(R*C),X,1,Xt,1);
        if (iscolmajor)
        {
            // This has ~identical speed, but code is longer
            // for (size_t r=0; r<R; ++r, Xt-=R*C-1)
            // {
            //     for (size_t c=0; c<C; ++c, Xt+=R) { *X++ = *Xt; }
            // }
            for (size_t r=0; r<R; ++r, X+=C, ++Xt) { cblas_scopy((int)(C),Xt,(int)R,X,1); }
            Xt -= R;
        }
        else
        {
            for (size_t c=0; c<C; ++c, X+=R, ++Xt) { cblas_scopy((int)(R),Xt,(int)C,X,1); }
            Xt -= C;
        }
        free(Xt);
    }

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int transpose_inplace_d (double *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R==1 || C==1) {}
    else if (R==C)
    {
        for (size_t c=0, n=1; c<C-1; ++c, n+=R+1)
        {
            cblas_dswap((int)(R-c-1),&X[n],1,&X[(c+1)*R+c],(int)C);
        }
    }
    else
    {
        double *Xt;
        if (!(Xt=(double *)malloc(R*C*sizeof(double)))) { fprintf(stderr,"error in transpose_inplace_d: problem with malloc. "); perror("malloc"); return 1; }
        cblas_dcopy((int)(R*C),X,1,Xt,1);
        if (iscolmajor)
        {
            for (size_t r=0; r<R; ++r, X+=C, ++Xt) { cblas_dcopy((int)(C),Xt,(int)R,X,1); }
            Xt -= R;
        }
        else
        {
            for (size_t c=0; c<C; ++c, X+=R, ++Xt) { cblas_dcopy((int)(R),Xt,(int)C,X,1); }
            Xt -= C;
        }
        free(Xt);
    }

    return 0;
}


int transpose_inplace_c (float *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R==1 || C==1) {}
    else if (R==C)
    {
        for (size_t c=0, n=2; c<C-1; ++c, n+=2*R+2)
        {
            cblas_cswap((int)(R-c-1),&X[n],1,&X[2*((c+1)*R+c)],(int)C);
        }
    }
    else
    {
        float *Xt;
        if (!(Xt=(float *)malloc(2*R*C*sizeof(float)))) { fprintf(stderr,"error in transpose_inplace_c: problem with malloc. "); perror("malloc"); return 1; }
        cblas_ccopy((int)(R*C),X,1,Xt,1);
        if (iscolmajor)
        {
            for (size_t r=0; r<R; ++r, X+=2*C, Xt+=2) { cblas_ccopy((int)(C),Xt,(int)R,X,1); }
            Xt -= 2*R;
        }
        else
        {
            for (size_t c=0; c<C; ++c, X+=2*R, Xt+=2) { cblas_ccopy((int)(R),Xt,(int)C,X,1); }
            Xt -= 2*C;
        }
        free(Xt);
    }

    return 0;
}


int transpose_inplace_z (double *X, const size_t R, const size_t C, const char iscolmajor)
{
    if (R==1 || C==1) {}
    else if (R==C)
    {
        for (size_t c=0, n=2; c<C-1; ++c, n+=2*R+2)
        {
            cblas_zswap((int)(R-c-1),&X[n],1,&X[2*((c+1)*R+c)],(int)C);
        }
    }
    else
    {
        double *Xt;
        if (!(Xt=(double *)malloc(2*R*C*sizeof(double)))) { fprintf(stderr,"error in transpose_inplace_z: problem with malloc. "); perror("malloc"); return 1; }
        cblas_zcopy((int)(R*C),X,1,Xt,1);
        if (iscolmajor)
        {
            for (size_t r=0; r<R; ++r, X+=2*C, Xt+=2) { cblas_zcopy((int)(C),Xt,(int)R,X,1); }
            Xt -= 2*R;
        }
        else
        {
            for (size_t c=0; c<C; ++c, X+=2*R, Xt+=2) { cblas_zcopy((int)(R),Xt,(int)C,X,1); }
            Xt -= 2*C;
        }
        free(Xt);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
