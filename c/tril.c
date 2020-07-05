//Zeros all elements above the kth diagonal of matrix X.
//This has in-place and not-in-place versions.

//LAPACKE_?lacpy was surprisingly a bit slower even for square, k==0 case.

#include <stdio.h>
#include <cblas.h>
#include <lapacke.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int tril_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int tril_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int tril_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int tril_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);

int tril_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int tril_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int tril_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);
int tril_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k);


int tril_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_s: k must be in [1-R C-1]\n"); return 1; }

    const float z = 0.0f;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    // if (k==0 && R==C && S*H==1)
    // {
    //     const CBLAS_ORDER Ord = (iscolmajor) ? CblasColMajor : CblasRowMajor;
    //     const int lda = (iscolmajor) ? R : C;
    //     if (LAPACKE_slacpy (Ord,'L',(int)R,(int)C,X,lda,Y,lda))
    //     { fprintf(stderr,"error in tril_s: problem with LAPACKE function\n"); return 1; }
    // }
    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0;                           //number of all-X cols
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                if (CX>0) { cblas_scopy((int)(R*CX),&X[n],1,&Y[n],1); n += R*CX; }
                for (size_t c=CX; c<C-C0; c++)
                {
                    cblas_scopy((int)c-k,&z,0,&Y[n],1); n += (size_t)((int)c-k);
                    cblas_scopy((int)R-(int)c+k,&X[n],1,&Y[n],1); n += (size_t)((int)R-(int)c+k);
                }
                if (C0>0) { cblas_scopy((int)(R*C0),&z,0,&Y[n],1); n += R*C0; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0;    //number of all-X rows
        if (R0>0) { cblas_scopy((int)(R0*C*SH),&z,0,Y,1); }
        for (size_t r=R0, n=R0*C*SH; r<R-RX; r++)
        {
            cblas_scopy((int)SH*((int)r+k+1),&X[n],1,&Y[n],1); n += SH*(size_t)((int)r+k+1);
            cblas_scopy((int)SH*((int)C-(int)r-k-1),&z,0,&Y[n],1); n += SH*(size_t)((int)C-(int)r-k-1);
        }
        if (RX>0) { cblas_scopy((int)(RX*C*SH),&X[C*SH*(R-RX)],1,&Y[C*SH*(R-RX)],1); }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int tril_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_d: k must be in [1-R C-1]\n"); return 1; }

    const double z = 0.0;

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0;                           //number of all-X cols
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                if (CX>0) { cblas_dcopy((int)(R*CX),&X[n],1,&Y[n],1); n += R*CX; }
                for (size_t c=CX; c<C-C0; c++)
                {
                    cblas_dcopy((int)c-k,&z,0,&Y[n],1); n += (size_t)((int)c-k);
                    cblas_dcopy((int)R-(int)c+k,&X[n],1,&Y[n],1); n += (size_t)((int)R-(int)c+k);
                }
                if (C0>0) { cblas_dcopy((int)(R*C0),&z,0,&Y[n],1); n += R*C0; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0;    //number of all-X rows
        if (R0>0) { cblas_dcopy((int)(R0*C*SH),&z,0,Y,1); }
        for (size_t r=R0, n=R0*C*SH; r<R-RX; r++)
        {
            cblas_dcopy((int)SH*((int)r+k+1),&X[n],1,&Y[n],1); n += SH*(size_t)((int)r+k+1);
            cblas_dcopy((int)SH*((int)C-(int)r-k-1),&z,0,&Y[n],1); n += SH*(size_t)((int)C-(int)r-k-1);
        }
        if (RX>0) { cblas_dcopy((int)(RX*C*SH),&X[C*SH*(R-RX)],1,&Y[C*SH*(R-RX)],1); }
    }

    return 0;
}


int tril_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_c: k must be in [1-R C-1]\n"); return 1; }

    const float z[2] = {0.0f,0.0f};

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0;                           //number of all-X cols
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                if (CX>0) { cblas_ccopy((int)(R*CX),&X[n],1,&Y[n],1); n += 2*R*CX; }
                for (size_t c=CX; c<C-C0; c++)
                {
                    cblas_ccopy((int)c-k,z,0,&Y[n],1); n += 2*(size_t)((int)c-k);
                    cblas_ccopy((int)R-(int)c+k,&X[n],1,&Y[n],1); n += 2*(size_t)((int)R-(int)c+k);
                }
                if (C0>0) { cblas_ccopy((int)(R*C0),z,0,&Y[n],1); n += 2*R*C0; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0;    //number of all-X rows
        if (R0>0) { cblas_ccopy((int)(R0*C*SH),z,0,Y,1); }
        for (size_t r=R0, n=2*R0*C*SH; r<R-RX; r++)
        {
            cblas_ccopy((int)SH*((int)r+k+1),&X[n],1,&Y[n],1); n += 2*SH*(size_t)((int)r+k+1);
            cblas_ccopy((int)SH*((int)C-(int)r-k-1),z,0,&Y[n],1); n += 2*SH*(size_t)((int)C-(int)r-k-1);
        }
        if (RX>0) { cblas_ccopy((int)(RX*C*SH),&X[2*C*SH*(R-RX)],1,&Y[2*C*SH*(R-RX)],1); }
    }

    return 0;
}


int tril_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_z: k must be in [1-R C-1]\n"); return 1; }

    const double z[2] = {0.0,0.0};

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0;                           //number of all-X cols
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                if (CX>0) { cblas_zcopy((int)(R*CX),&X[n],1,&Y[n],1); n += 2*R*CX; }
                for (size_t c=CX; c<C-C0; c++)
                {
                    cblas_zcopy((int)c-k,z,0,&Y[n],1); n += 2*(size_t)((int)c-k);
                    cblas_zcopy((int)R-(int)c+k,&X[n],1,&Y[n],1); n += 2*(size_t)((int)R-(int)c+k);
                }
                if (C0>0) { cblas_zcopy((int)(R*C0),z,0,&Y[n],1); n += 2*R*C0; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0;    //number of all-X rows
        if (R0>0) { cblas_zcopy((int)(R0*C*SH),z,0,Y,1); }
        for (size_t r=R0, n=2*R0*C*SH; r<R-RX; r++)
        {
            cblas_zcopy((int)SH*((int)r+k+1),&X[n],1,&Y[n],1); n += 2*SH*(size_t)((int)r+k+1);
            cblas_zcopy((int)SH*((int)C-(int)r-k-1),z,0,&Y[n],1); n += 2*SH*(size_t)((int)C-(int)r-k-1);
        }
        if (RX>0) { cblas_zcopy((int)(RX*C*SH),&X[2*C*SH*(R-RX)],1,&Y[2*C*SH*(R-RX)],1); }
    }

    return 0;
}


int tril_inplace_s (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_inplace_s: k must be in [1-R C-1]\n"); return 1; }

    const float z = 0.0f;
    //struct timespec tic, toc;
    //clock_gettime(CLOCK_REALTIME,&tic);

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0;                           //number of all-X cols
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                if (CX>0) { n += R*CX; }
                for (size_t c=CX; c<C-C0; c++)
                {
                    cblas_scopy((int)c-k,&z,0,&X[n],1);
                    n += R;
                }
                if (C0>0) { cblas_scopy((int)(R*C0),&z,0,&X[n],1); n += R*C0; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0;    //number of all-X rows
        if (R0>0) { cblas_scopy((int)(R0*C*SH),&z,0,X,1); }
        for (size_t r=R0, n=R0*C*SH; r<R-RX; r++)
        {
            n += SH*(size_t)((int)r+k+1);
            cblas_scopy((int)SH*((int)C-(int)r-k-1),&z,0,&X[n],1);
            n += SH;
        }
    }
    //clock_gettime(CLOCK_REALTIME,&toc);
    //fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int tril_inplace_d (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_inplace_d: k must be in [1-R C-1]\n"); return 1; }

    const double z = 0.0;

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0;                           //number of all-X cols
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                if (CX>0) { n += R*CX; }
                for (size_t c=CX; c<C-C0; c++)
                {
                    cblas_dcopy((int)c-k,&z,0,&X[n],1);
                    n += R;
                }
                if (C0>0) { cblas_dcopy((int)(R*C0),&z,0,&X[n],1); n += R*C0; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0;    //number of all-X rows
        if (R0>0) { cblas_dcopy((int)(R0*C*SH),&z,0,X,1); }
        for (size_t r=R0, n=R0*C*SH; r<R-RX; r++)
        {
            n += SH*(size_t)((int)r+k+1);
            cblas_dcopy((int)SH*((int)C-(int)r-k-1),&z,0,&X[n],1);
            n += SH*(size_t)((int)C-(int)r-k-1);
        }
    }

    return 0;
}


int tril_inplace_c (float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_inplace_c: k must be in [1-R C-1]\n"); return 1; }

    const float z[2] = {0.0f,0.0f};

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0;                           //number of all-X cols
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                if (CX>0) { n += 2*R*CX; }
                for (size_t c=CX; c<C-C0; c++)
                {
                    cblas_ccopy((int)c-k,z,0,&X[n],1);
                    n += 2*R;
                }
                if (C0>0) { cblas_ccopy((int)(R*C0),z,0,&X[n],1); n += 2*R*C0; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0;    //number of all-X rows
        if (R0>0) { cblas_ccopy((int)(R0*C*SH),z,0,X,1); }
        for (size_t r=R0, n=2*R0*C*SH; r<R-RX; r++)
        {
            n += 2*SH*(size_t)((int)r+k+1);
            cblas_ccopy((int)SH*((int)C-(int)r-k-1),z,0,&X[n],1);
            n += 2*SH*(size_t)((int)C-(int)r-k-1);
        }
    }

    return 0;
}


int tril_inplace_z (double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const int k)
{
    if (k<=-(int)R || k>=(int)C) { fprintf(stderr,"error in tril_inplace_z: k must be in [1-R C-1]\n"); return 1; }

    const double z[2] = {0.0,0.0};

    if (iscolmajor)
    {
        const size_t C0 = (k<(int)C-(int)R) ? (size_t)((int)C-(int)R-k) : 0;    //number of all-0 cols
        const size_t CX = (k>-1) ? (size_t)(k+1) : 0;                           //number of all-X cols
        for (size_t h=0, n=0; h<H; h++)
        {
            for (size_t s=0; s<S; s++)
            {
                if (CX>0) { n += 2*R*CX; }
                for (size_t c=CX; c<C-C0; c++)
                {
                    cblas_zcopy((int)c-k,z,0,&X[n],1);
                    n += 2*R;
                }
                if (C0>0) { cblas_zcopy((int)(R*C0),z,0,&X[n],1); n += 2*R*C0; }
            }
        }
    }
    else
    {
        const size_t SH = S*H;
        const size_t R0 = (k<0) ? (size_t)(-k) : 0;                                 //number of all-0 rows
        const size_t RX = (k>(int)C-(int)R-1) ? (size_t)((int)R-(int)C+k+1) : 0;    //number of all-X rows
        if (R0>0) { cblas_zcopy((int)(R0*C*SH),z,0,X,1); }
        for (size_t r=R0, n=2*R0*C*SH; r<R-RX; r++)
        {
            n += 2*SH*(size_t)((int)r+k+1);
            cblas_zcopy((int)SH*((int)C-(int)r-k-1),z,0,&X[n],1);
            n += 2*SH*(size_t)((int)C-(int)r-k-1);
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
